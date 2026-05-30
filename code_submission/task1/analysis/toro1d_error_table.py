#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

VARS = ["rho", "u", "p", "eint"]


def read_csv(path):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        rows = list(r)
    out = {"x": []}
    for v in VARS:
        out[v] = []
    for row in rows:
        out["x"].append(float(row["x"]))
        for v in VARS:
            out[v].append(float(row[v]))
    return out


def cell_widths(xs):
    n = len(xs)
    if n <= 1:
        return [1.0]
    dx = [0.0] * n
    dx[0] = xs[1] - xs[0]
    dx[-1] = xs[-1] - xs[-2]
    for i in range(1, n - 1):
        dx[i] = 0.5 * (xs[i + 1] - xs[i - 1])
    return dx


def l1_error(a, b, w):
    num = 0.0
    den = 0.0
    for ai, bi, wi in zip(a, b, w):
        num += abs(ai - bi) * wi
        den += wi
    return num / den if den > 0.0 else float("nan")


def rate(e_coarse, e_fine):
    if e_coarse <= 0.0 or e_fine <= 0.0:
        return float("nan")
    return math.log(e_coarse / e_fine, 2.0)


def parse_csv_list(text, cast=int):
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser(description="Build Toro L1 error and convergence tables")
    ap.add_argument("--results-dir", default=".")
    ap.add_argument("--modes", default="uniform,amr")
    ap.add_argument("--tests", default="1,2,3,4,5")
    ap.add_argument("--res", default="100,200,400", help="effective resolutions")
    ap.add_argument("--out-csv", default="toro1d_errors.csv")
    ap.add_argument("--out-md", default="toro1d_errors.md")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    modes = [x.strip() for x in args.modes.split(",") if x.strip()]
    tests = parse_csv_list(args.tests, int)
    res = parse_csv_list(args.res, int)

    rows = []
    rates = {}

    for mode in modes:
        rates[mode] = {}
        for test in tests:
            rates[mode][test] = {}
            data = {}
            for n_eff in res:
                dpath = results_dir / f"{mode}_t{test}_n{n_eff}_derived.csv"
                epath = results_dir / f"{mode}_t{test}_n{n_eff}_exact.csv"
                if not dpath.exists() or not epath.exists():
                    raise SystemExit(f"Missing CSV pair: {dpath.name} / {epath.name}")
                d = read_csv(dpath)
                e = read_csv(epath)
                if len(d["x"]) != len(e["x"]):
                    raise SystemExit(f"Length mismatch for {dpath.name} vs {epath.name}")
                w = cell_widths(d["x"])
                errs = {v: l1_error(d[v], e[v], w) for v in VARS}
                data[n_eff] = errs
                row = {
                    "mode": mode,
                    "test": test,
                    "effective_n": n_eff,
                }
                row.update(errs)
                rows.append(row)

            rtab = {}
            for i in range(1, len(res)):
                n0 = res[i - 1]
                n1 = res[i]
                key = f"{n0}->{n1}"
                rtab[key] = {}
                for v in VARS:
                    rtab[key][v] = rate(data[n0][v], data[n1][v])
            rates[mode][test] = rtab

    out_csv = results_dir / args.out_csv
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "mode", "test", "effective_n", "rho", "u", "p", "eint",
            "order_from_prev", "ord_rho", "ord_u", "ord_p", "ord_eint"
        ])
        for row in sorted(rows, key=lambda r: (r["mode"], r["test"], r["effective_n"])):
            n_eff = row["effective_n"]
            idx = res.index(n_eff)
            if idx == 0:
                key = ""
                ords = ["", "", "", ""]
            else:
                key = f"{res[idx - 1]}->{res[idx]}"
                rr = rates[row["mode"]][row["test"]][key]
                ords = [rr[v] for v in VARS]
            w.writerow([
                row["mode"], row["test"], row["effective_n"],
                row["rho"], row["u"], row["p"], row["eint"],
                key, ords[0], ords[1], ords[2], ords[3],
            ])

    out_md = results_dir / args.out_md
    with open(out_md, "w") as f:
        for mode in modes:
            f.write(f"## Mode: {mode}\n\n")
            f.write("| Test | N | L1(rho) | L1(u) | L1(p) | L1(eint) | ord_rho | ord_u | ord_p | ord_eint |\n")
            f.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for test in tests:
                for idx, n_eff in enumerate(res):
                    row = next(r for r in rows if r["mode"] == mode and r["test"] == test and r["effective_n"] == n_eff)
                    if idx == 0:
                        ords = ["-", "-", "-", "-"]
                    else:
                        key = f"{res[idx - 1]}->{res[idx]}"
                        rr = rates[mode][test][key]
                        ords = [f"{rr[v]:.3f}" for v in VARS]
                    f.write(
                        f"| {test} | {n_eff} | {row['rho']:.6e} | {row['u']:.6e} | {row['p']:.6e} | {row['eint']:.6e} | "
                        f"{ords[0]} | {ords[1]} | {ords[2]} | {ords[3]} |\n"
                    )
                f.write("\n")

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
