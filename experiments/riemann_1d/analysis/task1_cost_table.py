#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


def parse_int_list(text):
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def parse_log_metrics(path):
    runtime = None
    total_updates = 0
    updates_by_level = {}

    rt_re = re.compile(r"Run time =\s*([0-9.eE+-]+)")
    upd_re = re.compile(r"\[Level\s+(\d+)\s+step\s+\d+\]\s+Advanced\s+(\d+)\s+cells")

    with open(path) as f:
        for line in f:
            mrt = rt_re.search(line)
            if mrt:
                runtime = float(mrt.group(1))
            mup = upd_re.search(line)
            if mup:
                lev = int(mup.group(1))
                cnt = int(mup.group(2))
                total_updates += cnt
                updates_by_level[lev] = updates_by_level.get(lev, 0) + cnt

    if runtime is None:
        raise RuntimeError(f"Could not parse runtime from {path}")

    return runtime, total_updates, updates_by_level


def safe_ratio(num, den):
    if den == 0:
        return float("nan")
    return num / den


def main():
    ap = argparse.ArgumentParser(description="Create Task-1 AMR vs uniform cost tables from existing logs.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--tests", default="1,2,3,4,5")
    ap.add_argument("--res", default="100,200,400")
    ap.add_argument("--out-csv", default="task1_cost_metrics.csv")
    ap.add_argument("--out-md", default="task1_cost_metrics.md")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    tests = parse_int_list(args.tests)
    res = parse_int_list(args.res)

    rows = []
    for t in tests:
        for n in res:
            uniform_log = results_dir / f"uniform_t{t}_n{n}.log"
            amr_log = results_dir / f"amr_t{t}_n{n}.log"
            if not uniform_log.exists() or not amr_log.exists():
                raise FileNotFoundError(f"Missing expected logs for test={t}, n={n}")

            u_rt, u_upd, u_lvl = parse_log_metrics(uniform_log)
            a_rt, a_upd, a_lvl = parse_log_metrics(amr_log)

            row = {
                "test": t,
                "effective_n": n,
                "runtime_uniform_s": u_rt,
                "runtime_amr_s": a_rt,
                "runtime_ratio_amr_over_uniform": safe_ratio(a_rt, u_rt),
                "updates_uniform_total": u_upd,
                "updates_amr_total": a_upd,
                "updates_ratio_amr_over_uniform": safe_ratio(a_upd, u_upd),
                "amr_updates_level0": a_lvl.get(0, 0),
                "amr_updates_level1": a_lvl.get(1, 0),
                "amr_updates_level2": a_lvl.get(2, 0),
            }
            rows.append(row)

    out_csv = results_dir / args.out_csv
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "test", "effective_n",
            "runtime_uniform_s", "runtime_amr_s", "runtime_ratio_amr_over_uniform",
            "updates_uniform_total", "updates_amr_total", "updates_ratio_amr_over_uniform",
            "amr_updates_level0", "amr_updates_level1", "amr_updates_level2",
        ])
        for r in rows:
            w.writerow([
                r["test"], r["effective_n"],
                r["runtime_uniform_s"], r["runtime_amr_s"], r["runtime_ratio_amr_over_uniform"],
                r["updates_uniform_total"], r["updates_amr_total"], r["updates_ratio_amr_over_uniform"],
                r["amr_updates_level0"], r["amr_updates_level1"], r["amr_updates_level2"],
            ])

    out_md = results_dir / args.out_md
    with open(out_md, "w") as f:
        f.write("# Task 1 Cost Metrics (from existing logs)\n\n")
        f.write(
            "| Test | N_eff | Runtime uniform (s) | Runtime AMR (s) | AMR/Uniform runtime | "
            "Updated cells uniform | Updated cells AMR | AMR/Uniform updates |\n"
        )
        f.write("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in rows:
            f.write(
                f"| {r['test']} | {r['effective_n']} | {r['runtime_uniform_s']:.6f} | "
                f"{r['runtime_amr_s']:.6f} | {r['runtime_ratio_amr_over_uniform']:.3f} | "
                f"{r['updates_uniform_total']} | {r['updates_amr_total']} | "
                f"{r['updates_ratio_amr_over_uniform']:.3f} |\n"
            )

        f.write("\n")
        f.write("## Resolution-wise mean ratios\n\n")
        for n_focus in res:
            focus = [r for r in rows if r["effective_n"] == n_focus]
            if not focus:
                continue
            avg_rt = sum(r["runtime_ratio_amr_over_uniform"] for r in focus) / len(focus)
            avg_upd = sum(r["updates_ratio_amr_over_uniform"] for r in focus) / len(focus)
            f.write(f"- N_eff={n_focus}: mean runtime ratio (AMR/uniform) = **{avg_rt:.3f}**, ")
            f.write(f"mean updated-cell ratio (AMR/uniform) = **{avg_upd:.3f}**\n")
        f.write("\n")
        f.write("- Interpretation: updated-cell count measures algorithmic work; runtime additionally includes AMR bookkeeping overhead.\n")

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
