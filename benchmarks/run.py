#!/usr/bin/env python3
"""Repeat isolated FERRE runs and compare ID-keyed OPF numerical output."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import statistics
import subprocess
import tempfile
import time


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def control(text, threads):
    # Deliberately require a single simple namelist, avoiding ambiguous edits.
    text = '\n'.join(line.split('!')[0] for line in text.splitlines())
    if len(re.findall(r'&lista\b', text, re.I)) != 1 or not text.rstrip().endswith('/'):
        raise ValueError('Expected one &LISTA namelist ending in /')
    text = re.sub(r'\bnthreads\s*=\s*\d+\s*,?', '', text, flags=re.I)
    return text.rstrip()[:-1] + f'\n NTHREADS = {threads}\n /\n'


def output_paths(text):
    paths = {}
    for key in ('opfile', 'offile', 'sffile'):
        m = re.search(r'\b' + key + r'\s*=\s*[\'\"]([^\'\"]+)[\'\"]', text, re.I)
        if m:
            p = Path(m.group(1).strip())
            if p.is_absolute() or '..' in p.parts:
                raise ValueError(f'{key} must stay inside the case directory')
            paths[key] = p
    if 'opfile' not in paths:
        raise ValueError('Set an explicit OPFILE for numerical validation')
    return paths


def read_opf(path):
    rows = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] in rows:
            raise ValueError('Duplicate OPF object ID: ' + fields[0])
        values = [float(v.replace('D', 'E').replace('d', 'e')) for v in fields[1:]]
        if not values or not all(math.isfinite(v) for v in values):
            raise ValueError('Invalid/nonfinite OPF row: ' + fields[0])
        rows[fields[0]] = values
    if not rows:
        raise ValueError('Empty OPF output')
    return rows


def compare(reference, candidate, atol, rtol):
    if reference.keys() != candidate.keys():
        return {'passed': False, 'reason': 'Object IDs differ'}
    maxdiff = []
    passed = True
    for key, a in reference.items():
        b = candidate[key]
        if len(a) != len(b):
            return {'passed': False, 'reason': 'Column counts differ'}
        if not maxdiff:
            maxdiff = [0.0] * len(a)
        for i, (x, y) in enumerate(zip(a, b)):
            maxdiff[i] = max(maxdiff[i], abs(x-y))
            passed &= math.isclose(x, y, abs_tol=atol, rel_tol=rtol)
    return {'passed': bool(passed), 'max_abs_difference_by_column': maxdiff}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--case', required=True, type=Path, help='Directory copied into an isolated workspace')
    p.add_argument('--nml', default='input.nml_all')
    p.add_argument('--threads', nargs='+', type=int, default=[1, 2, 4, 8])
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--warmups', type=int, default=1)
    p.add_argument('--timeout', type=float, default=1800)
    p.add_argument('--out', required=True, type=Path, help='New directory, outside case directory')
    p.add_argument('--baseline', type=Path, help='Previous results.json from the same case')
    p.add_argument('--atol', type=float, default=1e-5)
    p.add_argument('--rtol', type=float, default=1e-5)
    p.add_argument('--build-label', default='', help='Compiler/version/flags or production executable description')
    a = p.parse_args()
    if min(a.threads) < 1 or a.repeats < 1 or a.warmups < 0:
        p.error('Use positive threads/repeats and nonnegative warmups')
    exe, case, out = a.exe.resolve(), a.case.resolve(), a.out.resolve()
    if out == case or case in out.parents:
        p.error('--out must be outside --case')
    source = (case / a.nml).read_text()
    outputs = output_paths(source)
    # Hash all staged inputs to detect changed data/configurations across builds.
    inputs = {str(f.relative_to(case)): digest(f) for f in sorted(case.rglob('*'))
              if f.is_file() and f.relative_to(case) not in outputs.values()}
    report = {'schema_version': 1, 'executable_sha256': digest(exe),
              'build_label': a.build_label, 'platform': platform.platform(),
              'cpu_count': os.cpu_count(), 'inputs_sha256': inputs,
              'nml': a.nml, 'runs': [], 'summary': [],
              'timing_scope': 'Whole process: grid load, fitting, uncertainty estimation, output and sorting'}
    baseline = json.loads(a.baseline.read_text()) if a.baseline else None
    if baseline and (baseline['inputs_sha256'] != inputs or baseline['nml'] != a.nml):
        p.error('Baseline inputs/configuration differ')
    reference = baseline['reference_opf'] if baseline else None
    out.mkdir(parents=True, exist_ok=False)
    failed = False
    for threads in a.threads:
        times = []
        for iteration in range(a.warmups + a.repeats):
            label = f't{threads}-r{iteration}'
            with tempfile.TemporaryDirectory(prefix='ferre-bench-') as tmp:
                work = Path(tmp) / 'case'
                shutil.copytree(case, work)
                for path in outputs.values():
                    (work / path).unlink(missing_ok=True)
                (work / 'benchmark.nml').write_text(control(source, threads))
                env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC='FALSE',
                           OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
                start = time.perf_counter()
                with (out / (label + '.log')).open('w') as log:
                    try:
                        result = subprocess.run([str(exe), 'benchmark.nml'], cwd=work,
                                                env=env, stdout=log, stderr=subprocess.STDOUT,
                                                timeout=a.timeout)
                    except subprocess.TimeoutExpired:
                        raise SystemExit(f'{label} timed out; see {out}')
                elapsed = time.perf_counter() - start
                if result.returncode:
                    raise SystemExit(f'{label} failed ({result.returncode}); see {out}')
                # Fortran STOP can return zero on errors: require valid numerical output too.
                rows = read_opf(work / outputs['opfile'])
                if reference is None:
                    reference = rows
                check = compare(reference, rows, a.atol, a.rtol)
                failed |= not check['passed']
                for key, path in outputs.items():
                    if (work / path).is_file():
                        shutil.copy2(work / path, out / (label + '.' + key))
                warmup = iteration < a.warmups
                report['runs'].append({'threads': threads, 'warmup': warmup,
                                       'seconds': elapsed, 'objects': len(rows), 'comparison': check})
                if not warmup:
                    times.append(elapsed)
                print(f'{label}: {elapsed:.3f}s; agreement={check["passed"]}', flush=True)
        median = statistics.median(times)
        report['summary'].append({'threads': threads, 'median_seconds': median,
                                  'min_seconds': min(times), 'max_seconds': max(times),
                                  'objects_per_hour': len(rows) * 3600 / median})
    report['reference_opf'] = reference
    report['agreement_passed'] = not failed
    (out / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['summary'], indent=2))
    raise SystemExit(1 if failed else 0)


if __name__ == '__main__':
    main()
