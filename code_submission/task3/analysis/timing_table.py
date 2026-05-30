#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def to_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def to_int(x):
    try:
        return int(x)
    except Exception:
        return -1


def main():
    ap = argparse.ArgumentParser(description="Summarize quadrant timing study CSV")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out-md", default="timing_summary.md")
    args = ap.parse_args()

    rows = read_rows(args.csv)

    # Normalize numeric fields for sorting/metrics.
    for r in rows:
        r["np_i"] = to_int(r.get("np", "1"))
        r["eff_i"] = to_int(r.get("effective_n", "-1"))
        r["rt"] = to_float(r.get("runtime_s", "nan"))

    out_md = Path(args.out_md)
    with open(out_md, "w") as f:
        f.write("## Raw Timing Table\n\n")
        f.write("| case_id | group | eff_N | max_level | np | runtime_s |\n")
        f.write("|---|---|---:|---:|---:|---:|\n")
        for r in sorted(rows, key=lambda x: (x["group"], x["eff_i"], x["np_i"], x["case_id"])):
            f.write(
                f"| {r['case_id']} | {r['group']} | {r['eff_i']} | {r.get('max_level','')} | {r['np_i']} | {r['rt']:.6f} |\n"
            )

        # Uniform scaling tables by effective resolution.
        f.write("\n## Uniform MPI Scaling\n\n")
        uniform = [r for r in rows if r["group"] == "uniform_mpi"]
        by_eff = defaultdict(list)
        for r in uniform:
            by_eff[r["eff_i"]].append(r)

        for eff in sorted(by_eff.keys()):
            data = sorted(by_eff[eff], key=lambda x: x["np_i"])
            base = next((x for x in data if x["np_i"] == 1), None)
            if base is None or base["rt"] <= 0.0:
                continue
            f.write(f"### Uniform N={eff}\n\n")
            f.write("| np | runtime_s | speedup | parallel_efficiency |\n")
            f.write("|---:|---:|---:|---:|\n")
            for r in data:
                s = base["rt"] / r["rt"] if r["rt"] > 0.0 else float("nan")
                e = s / r["np_i"] if r["np_i"] > 0 else float("nan")
                f.write(f"| {r['np_i']} | {r['rt']:.6f} | {s:.3f} | {e:.3f} |\n")
            f.write("\n")

        # AMR high-resolution scaling.
        f.write("## AMR Highest-Resolution MPI Scaling\n\n")
        amr_hi = [r for r in rows if r["group"] == "amr_mpi_high"]
        amr_hi = sorted(amr_hi, key=lambda x: x["np_i"])
        base_hi = next((x for x in amr_hi if x["np_i"] == 1), None)
        if base_hi is not None and base_hi["rt"] > 0.0:
            f.write("| np | runtime_s | speedup | parallel_efficiency |\n")
            f.write("|---:|---:|---:|---:|\n")
            for r in amr_hi:
                s = base_hi["rt"] / r["rt"] if r["rt"] > 0.0 else float("nan")
                e = s / r["np_i"] if r["np_i"] > 0 else float("nan")
                f.write(f"| {r['np_i']} | {r['rt']:.6f} | {s:.3f} | {e:.3f} |\n")

        # Best overall speed-up requested in assignment:
        # compare highest uniform resolution at np=1 vs fastest AMR high (any np).
        f.write("\n## Best High-Resolution Speed-Up\n\n")
        uniform_hi_candidates = [r for r in rows if r["group"] in ("uniform_serial", "uniform_mpi")]
        if uniform_hi_candidates:
            max_eff = max(r["eff_i"] for r in uniform_hi_candidates)
            uni_hi_1 = next((r for r in rows if r["eff_i"] == max_eff and r["np_i"] == 1 and r["group"] in ("uniform_serial", "uniform_mpi")), None)
            amr_hi_any = [r for r in rows if r["group"] == "amr_mpi_high" and r["rt"] > 0.0]
            if uni_hi_1 is not None and amr_hi_any:
                best_amr = min(amr_hi_any, key=lambda x: x["rt"])
                sp = uni_hi_1["rt"] / best_amr["rt"]
                f.write(
                    f"Uniform high (N={max_eff}, np=1) runtime = {uni_hi_1['rt']:.6f} s\\\n\n"
                    f"Fastest AMR high runtime = {best_amr['rt']:.6f} s (np={best_amr['np_i']})\\\n\n"
                    f"Best speed-up = {sp:.3f}x\n"
                )

    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
