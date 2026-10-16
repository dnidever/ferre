import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('runner', Path(__file__).with_name('run.py'))
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)


class RunnerTests(unittest.TestCase):
    def test_order_and_tolerance(self):
        self.assertTrue(r.compare({'a': [1.], 'b': [2.]}, {'b': [2.], 'a': [1.]}, 0, 0)['passed'])
        self.assertFalse(r.compare({'a': [1.]}, {'a': [2.]}, 1e-5, 1e-5)['passed'])
        self.assertFalse(r.compare({'a': [1.]}, {'b': [1.]}, 0, 0)['passed'])

    def test_control_and_paths(self):
        s = "&LISTA\n NTHREADS=2,\n OPFILE='out.opf'\n /"
        self.assertEqual(r.control(s, 4).upper().count('NTHREADS'), 1)
        with self.assertRaises(ValueError):
            r.output_paths("OPFILE='../escape'")

    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case = root / 'case'
            case.mkdir()
            (case / 'input.nml_all').write_text("&LISTA\n OPFILE='fit.opf'\n /\n")
            exe = root / 'fake'
            exe.write_text('#!' + sys.executable + "\nfrom pathlib import Path\nPath('fit.opf').write_text('star 1.0 2.0\\n')\n")
            exe.chmod(0o755)
            cmd = [sys.executable, str(Path(r.__file__)), '--exe', str(exe), '--case', str(case),
                   '--threads', '1', '2', '--warmups', '0', '--repeats', '2', '--out', str(root / 'out')]
            subprocess.run(cmd, check=True, capture_output=True)
            self.assertFalse((case / 'fit.opf').exists())
            subprocess.run(cmd[:-1] + [str(root / 'out2'), '--baseline', str(root / 'out/results.json')],
                           check=True, capture_output=True)
            exe.write_text(exe.read_text().replace('1.0 2.0', '9.0 2.0'))
            result = subprocess.run(cmd[:-1] + [str(root / 'out3'), '--baseline', str(root / 'out/results.json')],
                                    capture_output=True)
            self.assertEqual(result.returncode, 1)


if __name__ == '__main__':
    unittest.main()
