#!/usr/bin/env python3
import argparse
from pathlib import Path

from postprocess_1d import (
    check_x_coords_against_geometry,
    infer_plotfile_from_slice,
    parse_inputs,
    parse_plotfile_header,
    read_slice,
    wave_alignment_from_rows,
    write_exact,
    run_fextract,
)


def main():
    parser = argparse.ArgumentParser(description="Write exact Riemann reference on the plotfile x-grid")
    parser.add_argument("--plotfile", type=str, default=None, help="Plotfile directory (pltXXXXX)")
    parser.add_argument("--slice", type=str, default=None, help="Optional precomputed slice file")
    parser.add_argument("--inputs", type=str, default="inputs", help="Inputs file")
    parser.add_argument("--out", type=str, required=True, help="Output exact CSV")
    parser.add_argument("--time", type=float, default=None, help="Override time for exact solution")
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    inputs_path = (base / args.inputs).resolve()
    prob = parse_inputs(inputs_path)

    plotfile_path = Path(args.plotfile) if args.plotfile else None
    if args.slice:
        slice_path = Path(args.slice)
        if plotfile_path is None:
            plotfile_path = infer_plotfile_from_slice(slice_path)
    else:
        if plotfile_path is None:
            raise SystemExit("Provide --plotfile or --slice")
        slice_path = Path(str(plotfile_path) + ".slice")
        run_fextract(plotfile_path, slice_path)

    t_slice, rows = read_slice(slice_path)
    if not rows:
        raise SystemExit(f"No rows in {slice_path}")

    geom = None
    if plotfile_path is not None and plotfile_path.is_dir():
        geom = parse_plotfile_header(plotfile_path)
        x_err = check_x_coords_against_geometry([r[0] for r in rows], geom)
        print("Geometry check:"
              f" prob_lo={geom['prob_lo']} prob_hi={geom['prob_hi']}"
              f" ncell={geom['ncell']} dx={geom['dx']}"
              f" max|x_slice-x_cellcenter|={x_err:.3e}")

    t = args.time
    if t is None:
        if geom is not None:
            t = geom["time"]
        elif t_slice is not None:
            t = t_slice
        else:
            t = prob.get("stop_time", None)
    if t is None:
        raise SystemExit("Could not determine time for exact solution")

    xs = [r[0] for r in rows]
    out_path = Path(args.out)
    waves = write_exact(out_path, xs, t, prob["test"], prob["gamma"], prob["x0"])

    print(f"Exact setup: t={t} gamma={prob['gamma']} x0={prob['x0']} test={prob['test']}")
    print("Predicted wave locations (name, speed, x):")
    for name, speed, xpos in waves:
        print(f"  {name:24s} speed={speed: .12e} x={xpos: .12e}")
    aligns = wave_alignment_from_rows(rows, waves)
    print("Nearest numerical rho-gradient edges to predicted waves:")
    for name, _speed, xpos, xnum, dxabs in aligns:
        print(f"  {name:24s} x_exact={xpos: .12e} x_num~={xnum: .12e} |dx|={dxabs: .12e}")

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
