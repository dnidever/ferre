"""Convert binary PCA grids between float32 and native uint16 storage.

NumPy and f90nml are needed only for conversion (``pip install ferre[quantize]``),
not for running the compiled FERRE reader. No source files are modified.
"""

from pathlib import Path
import hashlib
import json
import math
import shutil
import tempfile


def _dependencies():
    try:
        import numpy as np
        import f90nml
    except ImportError as error:
        raise ImportError(
            'Install the converter dependencies with pip install "ferre[quantize]"'
        ) from error
    return np, f90nml


def _root(stream):
    """Read precisely one namelist, leaving the rest of the header untouched."""
    _, f90nml = _dependencies()
    text = bytearray()
    quote = None
    comment = False
    while True:
        b = stream.read(1)
        if not b:
            raise ValueError("Unterminated SYNTH header")
        text.extend(b)
        c = b.decode("ascii")
        if comment:
            if c == "\n":
                comment = False
        elif quote:
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c == "!":
            comment = True
        elif c == "/":
            return dict(f90nml.reads(text.decode("ascii"))["synth"])


def _layout(header):
    np, _ = _dependencies()
    header = Path(header).resolve()
    if header.suffix.lower() != ".hdr":
        raise ValueError("Supply a binary PCA .hdr library")
    with header.open("rb") as stream:
        root = _root(stream)
    blocks = np.atleast_1d(root.get("npca", [])).astype(int)
    ncoeff = int(root["npix"])
    if not len(blocks) or np.any(blocks <= 0) or ncoeff < 1 or ncoeff % len(blocks):
        raise ValueError("Expected PCA coefficients, not uncompressed spectra")
    transposed = int(root.get("transposed", 0))
    if transposed not in (0, 1):
        raise ValueError("TRANSPOSED must be 0 or 1")
    nodepath = None
    if root.get("nodefile"):
        nodepath = Path(root["nodefile"])
        if not nodepath.is_absolute():
            nodepath = header.parent / nodepath
        records = 0
        with nodepath.open() as stream:
            for line in stream:
                row = line.split("#", 1)[0].split("!", 1)[0].replace(",", " ").split()
                if not row:
                    continue
                if len(row) != int(root["n_of_dim"]) or not all(
                    math.isfinite(float(v.replace("D", "E").replace("d", "e")))
                    for v in row
                ):
                    raise ValueError("Invalid NODEFILE coordinate row")
                records += 1
    else:
        counts = np.atleast_1d(root["n_p"])
        if len(counts) != int(root["n_of_dim"]) or any(
            int(n) != n or n < 1 for n in counts
        ):
            raise ValueError("Invalid regular-grid dimensions")
        records = math.prod(map(int, counts))
    if records < 1:
        raise ValueError("Grid contains no records")
    return header, root, nodepath, records, ncoeff, transposed


def _digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_records(stream, values, first, records, transposed):
    if not transposed:
        stream.write(values.tobytes(order="C"))
    else:
        for column in range(values.shape[1]):
            stream.seek(values.dtype.itemsize * (column * records + first))
            stream.write(values[:, column].tobytes())


def _convert(source, directory, batch_size, decode, progress):
    np, f90nml = _dependencies()
    if type(batch_size) is not int or batch_size < 1:
        raise ValueError("batch_size must be positive")
    source, root, nodes, records, ncoeff, transposed = _layout(source)
    mode = str(root.get("coeff_format", "float32")).lower()
    if mode != ("uint16" if decode else "float32"):
        raise ValueError("Unexpected source COEFF_FORMAT")
    destination = Path(directory).expanduser().resolve()
    if destination.exists():
        raise FileExistsError(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = source.with_suffix(".u16" if decode else ".unf")
    dtype = "<u2" if decode else "<f4"
    expected = records * ncoeff * np.dtype(dtype).itemsize
    if payload.stat().st_size != expected:
        raise ValueError("Source payload size mismatch")
    shape = (ncoeff, records) if transposed else (records, ncoeff)
    data = np.memmap(payload, dtype=dtype, mode="r", shape=shape)

    def batch(first):
        end = min(records, first + batch_size)
        return np.asarray(
            data[:, first:end].T if transposed else data[first:end], dtype=np.float64
        )

    if decode:
        tablepath = Path(root["quantfile"])
        if not tablepath.is_absolute():
            tablepath = source.parent / tablepath
        table = dict(f90nml.read(tablepath)["quantization"])
        offset = np.atleast_1d(table["offset"]).astype(float)
        scale = np.atleast_1d(table["scale"]).astype(float)
        if (
            table.get("version") != 1
            or table.get("ncoeff") != ncoeff
            or offset.shape != (ncoeff,)
            or scale.shape != (ncoeff,)
            or not np.isfinite(offset).all()
            or not np.isfinite(scale).all()
            or np.any(scale <= 0)
        ):
            raise ValueError("Invalid quantization table")
    else:
        offset = np.full(ncoeff, np.inf)
        maximum = np.full(ncoeff, -np.inf)
        for first in range(0, records, batch_size):
            values = batch(first)
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite source coefficients")
            offset = np.minimum(offset, values.min(axis=0))
            maximum = np.maximum(maximum, values.max(axis=0))
            if progress:
                progress("calibrate", min(records, first + batch_size), records)
        scale = np.where(maximum > offset, (maximum - offset) / 65535.0, 1.0)
    if (
        np.any(np.abs(offset) > np.finfo(np.float32).max)
        or not np.isfinite(offset + scale * 65535).all()
        or np.any(np.abs(offset + scale * 65535) > np.finfo(np.float32).max)
    ):
        raise ValueError("Quantization endpoints exceed float32 range")
    stage = Path(
        tempfile.mkdtemp(prefix=destination.name + ".partial-", dir=destination.parent)
    )
    try:
        newroot = dict(root)
        if nodes:
            newroot["nodefile"] = source.stem + ".nodes"
            shutil.copy2(nodes, stage / newroot["nodefile"])
        if decode:
            newroot.pop("coeff_format", None)
            newroot.pop("quantfile", None)
        else:
            newroot.update(coeff_format="uint16", quantfile=source.stem + ".quant")
            f90nml.Namelist(
                quantization=dict(
                    version=1,
                    ncoeff=ncoeff,
                    offset=offset.tolist(),
                    scale=scale.tolist(),
                )
            ).write(stage / newroot["quantfile"])
        with source.open("rb") as original, (stage / source.name).open("wb") as output:
            _root(original)
            output.write(str(f90nml.Namelist(synth=newroot)).encode("ascii"))
            shutil.copyfileobj(original, output)
        result = stage / (source.stem + (".unf" if decode else ".u16"))
        with result.open("wb") as stream:
            stream.truncate(records * ncoeff * (4 if decode else 2))
            for first in range(0, records, batch_size):
                values = batch(first)
                if decode:
                    values = (offset + values * scale).astype("<f4")
                else:
                    values = np.rint((values - offset) / scale)
                    if np.any((values < 0) | (values > 65535)):
                        raise ValueError(
                            "Quantization overflow; no clipping is allowed"
                        )
                    values = values.astype("<u2")
                if not np.isfinite(values).all():
                    raise ValueError("Nonfinite output coefficients")
                _write_records(stream, values, first, records, transposed)
                if progress:
                    progress("write", min(records, first + batch_size), records)
        files = {
            p.name: dict(bytes=p.stat().st_size, sha256=_digest(p))
            for p in stage.iterdir()
        }
        manifest = dict(
            storage="float32" if decode else "uint16",
            format_version=1,
            byte_order="little",
            records=records,
            ncoeff=ncoeff,
            transposed=transposed,
            source_header=str(source),
            source_payload_sha256=_digest(payload),
            files=files,
            decode="float64 offset + scale * unsigned code, then round to float32",
        )
        (stage / "quantization.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(destination)
    except BaseException:
        shutil.rmtree(stage)
        raise
    return destination / source.name


def quantize_grid(source, directory, *, batch_size=1024, progress=None):
    """Write a native uint16 PCA bundle, retaining coordinate and wavelength order.

    Supports regular or NODEFILE grids, MULTI wavelength headers and both
    record layouts. Source coefficients must be little-endian float32 in the
    companion .unf file. All coefficients must be finite. The destination must
    not exist. Calibration and writing are batched; no decoded whole grid is
    allocated. A callback receives (phase, completed_records, total_records).
    Requires NumPy and f90nml. Returns the new header Path.
    """
    return _convert(source, directory, batch_size, False, progress)


def dequantize_grid(source, directory, *, batch_size=1024, progress=None):
    """Expand a native uint16 bundle for older FERRE versions, in batches.

    Restores float32 coefficients and removes the quantization header fields.
    No interpolation or new PCA truncation is performed. The destination must
    not exist. Returns the compatible float32 header Path.
    """
    return _convert(source, directory, batch_size, True, progress)


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["quantize", "decode"])
    parser.add_argument("header", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--batch-size", type=int, default=1024)
    args = parser.parse_args()
    method = quantize_grid if args.mode == "quantize" else dequantize_grid
    print(method(args.header, args.directory, batch_size=args.batch_size))


if __name__ == "__main__":
    main()
