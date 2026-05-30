#!/usr/bin/env python3
import argparse
import csv
import math
import os
import platform
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def f64(x, default=math.nan):
    try:
        return float(x)
    except Exception:
        return default


def i64(x, default=-1):
    try:
        return int(x)
    except Exception:
        return default


def mean(vals):
    return sum(vals) / len(vals) if vals else math.nan


def std(vals):
    if len(vals) < 2:
        return 0.0
    return statistics.stdev(vals)


def ensure_dirs(results_dir: Path):
    (results_dir / "plots").mkdir(parents=True, exist_ok=True)
    (results_dir / "tables").mkdir(parents=True, exist_ok=True)


def save_fig(fig, base: Path):
    fig.savefig(base.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def aggregate_rows(rows, key_fields, value_field):
    out = defaultdict(list)
    for r in rows:
        key = tuple(r[k] for k in key_fields)
        out[key].append(f64(r[value_field]))
    return out


def compute_speedup_efficiency(mode_name, grouped_by_n_core):
    speedup_rows = []
    efficiency_rows = []
    for nval in sorted({k[0] for k in grouped_by_n_core.keys()}, key=lambda x: int(x)):
        by_core = {int(core): vals for (n, core), vals in grouped_by_n_core.items() if n == nval}
        if 1 not in by_core:
            continue
        t1 = mean(by_core[1])
        if not (t1 > 0):
            continue
        for core in sorted(by_core.keys()):
            t = mean(by_core[core])
            if not (t > 0):
                continue
            s = t1 / t
            e = s / core
            speedup_rows.append(
                {
                    "mode": mode_name,
                    "N_or_effectiveN": nval,
                    "cores": core,
                    "walltime_sec": t,
                    "speedup": s,
                }
            )
            efficiency_rows.append(
                {
                    "mode": mode_name,
                    "N_or_effectiveN": nval,
                    "cores": core,
                    "walltime_sec": t,
                    "efficiency": e,
                }
            )
    return speedup_rows, efficiency_rows


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def plot_time_vs_n(results_dir, uniform_rows):
    grouped = defaultdict(list)
    for r in uniform_rows:
        n = i64(r["N"])
        grouped[n].append(f64(r["walltime_sec"]))

    xs = sorted(grouped.keys())
    ys = [mean(grouped[x]) for x in xs]
    es = [std(grouped[x]) for x in xs]

    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.errorbar(xs, ys, yerr=es, marker="o", capsize=4, linewidth=2)
    ax.set_xlabel("Uniform resolution N (N x N)")
    ax.set_ylabel("Walltime (s)")
    ax.set_title("Task 3: Uniform No-AMR Serial Runtime vs Resolution")
    ax.grid(True, alpha=0.25)
    save_fig(fig, results_dir / "plots" / "time_vs_N")


def plot_time_vs_cores(results_dir, mpi_rows, amr_mpi_rows):
    grouped_uni = defaultdict(list)
    for r in mpi_rows:
        n = i64(r["N"])
        c = i64(r["cores"])
        grouped_uni[(n, c)].append(f64(r["walltime_sec"]))

    grouped_amr = defaultdict(list)
    for r in amr_mpi_rows:
        n = i64(r["effective_N"])
        c = i64(r["cores"])
        grouped_amr[(n, c)].append(f64(r["walltime_sec"]))

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    for n in sorted({k[0] for k in grouped_uni.keys()}):
        cores = sorted({k[1] for k in grouped_uni.keys() if k[0] == n})
        vals = [mean(grouped_uni[(n, c)]) for c in cores]
        ax.plot(cores, vals, marker="o", linewidth=2, label=f"Uniform N={n}")

    for n in sorted({k[0] for k in grouped_amr.keys()}):
        cores = sorted({k[1] for k in grouped_amr.keys() if k[0] == n})
        vals = [mean(grouped_amr[(n, c)]) for c in cores]
        ax.plot(cores, vals, marker="s", linestyle="--", linewidth=2, label=f"AMR eff N={n}")

    ax.set_xlabel("MPI cores")
    ax.set_ylabel("Walltime (s)")
    ax.set_title("Task 3: Runtime vs Cores")
    ax.set_xticks(sorted({i64(r["cores"]) for r in mpi_rows + amr_mpi_rows}))
    ax.grid(True, alpha=0.25)
    ax.legend()
    save_fig(fig, results_dir / "plots" / "time_vs_cores")


def plot_speedup_vs_cores(results_dir, speedup_rows):
    by_mode_n = defaultdict(list)
    for r in speedup_rows:
        by_mode_n[(r["mode"], str(r["N_or_effectiveN"]))].append(r)

    fig, ax = plt.subplots(figsize=(7.8, 5.0))
    max_core = 1
    for (mode, n), rows in sorted(by_mode_n.items(), key=lambda x: (x[0][0], int(x[0][1]))):
        rows = sorted(rows, key=lambda rr: i64(rr["cores"]))
        cores = [i64(rr["cores"]) for rr in rows]
        vals = [f64(rr["speedup"]) for rr in rows]
        max_core = max(max_core, max(cores))
        label = f"{mode} N={n}" if mode == "uniform_mpi" else f"{mode} eff N={n}"
        ax.plot(cores, vals, marker="o", linewidth=2, label=label)

    ideal_x = list(range(1, max_core + 1))
    ideal_y = ideal_x
    ax.plot(ideal_x, ideal_y, linestyle=":", color="black", label="Ideal linear")
    ax.set_xlabel("MPI cores")
    ax.set_ylabel("Speedup T(1)/T(p)")
    ax.set_title("Task 3: Speedup vs Cores")
    ax.grid(True, alpha=0.25)
    ax.legend()
    save_fig(fig, results_dir / "plots" / "speedup_vs_cores")


def plot_amr_vs_uniform(results_dir, uniform_rows, amr_match_rows):
    grouped_uni = defaultdict(list)
    for r in uniform_rows:
        grouped_uni[i64(r["N"])].append(f64(r["walltime_sec"]))

    grouped_amr = defaultdict(list)
    for r in amr_match_rows:
        eff = i64(r["effective_N"])
        grouped_amr[eff].append(f64(r["walltime_sec"]))

    targets = [1536, 3072]
    labels = [f"N={n}" for n in targets]
    uni_vals = [mean(grouped_uni[n]) if n in grouped_uni else math.nan for n in targets]
    uni_errs = [std(grouped_uni[n]) if n in grouped_uni else 0.0 for n in targets]
    amr_vals = [mean(grouped_amr[n]) if n in grouped_amr else math.nan for n in targets]
    amr_errs = [std(grouped_amr[n]) if n in grouped_amr else 0.0 for n in targets]

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    x = list(range(len(targets)))
    w = 0.35
    ax.bar(
        [v - w / 2 for v in x],
        uni_vals,
        yerr=uni_errs,
        capsize=4,
        width=w,
        label="Uniform serial",
        color="#4C78A8",
    )
    ax.bar(
        [v + w / 2 for v in x],
        amr_vals,
        yerr=amr_errs,
        capsize=4,
        width=w,
        label="AMR serial (base 768)",
        color="#F58518",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Walltime (s)")
    ax.set_title("Task 3: AMR Matching vs Uniform Runtime")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    save_fig(fig, results_dir / "plots" / "AMR_vs_uniform_time")


def plot_sensitivity(results_dir, sens_rows):
    rows = sorted(sens_rows, key=lambda r: r["case_id"])
    labels = [r["case_id"] for r in rows]
    runtimes = [f64(r["walltime_sec"]) for r in rows]
    lev1_cov = [f64(r.get("level1_coverage_pct", "nan")) for r in rows]

    fig, ax1 = plt.subplots(figsize=(8.6, 4.9))
    x = list(range(len(rows)))
    ax1.bar(x, runtimes, color="#54A24B", alpha=0.85)
    ax1.set_ylabel("Walltime (s)")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=25, ha="right")
    ax1.set_title("Task 3: Sensitivity Runtime and Level-1 Coverage")
    ax1.grid(True, axis="y", alpha=0.25)

    ax2 = ax1.twinx()
    ax2.plot(x, lev1_cov, color="#E45756", marker="o", linewidth=2)
    ax2.set_ylabel("Level-1 coverage (%)")

    save_fig(fig, results_dir / "plots" / "sensitivity_summary")


def write_summary_md(
    summary_path: Path,
    uniform_rows,
    mpi_rows,
    amr_match_rows,
    amr_mpi_rows,
    speedup_rows,
    efficiency_rows,
):
    grouped_uni = aggregate_rows(uniform_rows, ["N"], "walltime_sec")
    grouped_mpi = aggregate_rows(mpi_rows, ["N", "cores"], "walltime_sec")
    grouped_amr_mpi = aggregate_rows(amr_mpi_rows, ["effective_N", "cores"], "walltime_sec")

    lines = []
    lines.append("# Task 3 Quadrant Timing Summary")
    lines.append("")
    lines.append("## Fairness Precautions")
    lines.append("")
    lines.append(f"- Machine: `{platform.node()}`")
    lines.append(f"- Platform: `{platform.platform()}`")
    lines.append("- All timings recorded on the same laptop.")
    lines.append("- Quiet timing flags used (`adv.v=0`, `amr.v=0`, `amr.plot_int=1000000`).")
    lines.append("- MPI runs used `mpiexec -n p --bind-to none --map-by :OVERSUBSCRIBE`, matching the Task 3 timing scripts and writeup.")
    lines.append("- Laptop timings; scaling saturates due to hardware limits/thermals; all comparisons done on same machine.")
    lines.append("")

    lines.append("## Uniform Serial (No AMR)")
    lines.append("")
    lines.append("| N | repeats | mean_s | median_s | std_s |")
    lines.append("|---:|---:|---:|---:|---:|")
    for n in sorted({int(k[0]) for k in grouped_uni.keys()}):
        vals = grouped_uni[(str(n),)]
        lines.append(
            f"| {n} | {len(vals)} | {mean(vals):.6f} | {statistics.median(vals):.6f} | {std(vals):.6f} |"
        )
    lines.append("")

    lines.append("## Uniform MPI Speedup")
    lines.append("")
    lines.append("| N | cores | repeats | mean_s | median_s | std_s | speedup | efficiency |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
    uni_speed = [
        r for r in speedup_rows if r["mode"] == "uniform_mpi"
    ]
    uni_eff = {
        (r["N_or_effectiveN"], int(r["cores"])): r["efficiency"]
        for r in efficiency_rows
        if r["mode"] == "uniform_mpi"
    }
    for r in sorted(uni_speed, key=lambda rr: (int(rr["N_or_effectiveN"]), int(rr["cores"]))):
        n = int(r["N_or_effectiveN"])
        c = int(r["cores"])
        e = uni_eff[(r["N_or_effectiveN"], c)]
        vals = grouped_mpi.get((str(n), str(c)), [])
        rep = len(vals)
        med = statistics.median(vals) if vals else math.nan
        sd = std(vals)
        lines.append(
            f"| {n} | {c} | {rep} | {f64(r['walltime_sec']):.6f} | {med:.6f} | {sd:.6f} | {f64(r['speedup']):.3f} | {f64(e):.3f} |"
        )
    lines.append("")

    if amr_match_rows:
        lines.append("## AMR Match (Serial) Coverage")
        lines.append("")
        lines.append("| effective_N | repeats | mean_s | median_s | std_s | representative L0_cov_% | representative L1_cov_% | representative L2_cov_% |")
        lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
        by_eff = defaultdict(list)
        for r in amr_match_rows:
            by_eff[int(r["effective_N"])].append(r)
        for eff_n in sorted(by_eff.keys()):
            rows = by_eff[eff_n]
            vals = [f64(r["walltime_sec"]) for r in rows]
            cov_row = next(
                (
                    r
                    for r in rows
                    if not math.isnan(f64(r.get("level1_coverage_pct", "")))
                    or not math.isnan(f64(r.get("level2_coverage_pct", "")))
                ),
                rows[0],
            )
            lines.append(
                f"| {eff_n} | {len(vals)} | {mean(vals):.6f} | {statistics.median(vals):.6f} | {std(vals):.6f} | "
                f"{f64(cov_row.get('level0_coverage_pct', '')):.3f} | {f64(cov_row.get('level1_coverage_pct', '')):.3f} | {f64(cov_row.get('level2_coverage_pct', '')):.3f} |"
            )
        lines.append("")

    if grouped_amr_mpi:
        lines.append("## AMR MPI High-Resolution Speedup")
        lines.append("")
        lines.append("| effective_N | cores | repeats | mean_s | median_s | std_s | speedup | efficiency |")
        lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|")
        amr_speed = [r for r in speedup_rows if r["mode"] == "amr_mpi_high"]
        amr_eff = {
            (r["N_or_effectiveN"], int(r["cores"])): r["efficiency"]
            for r in efficiency_rows
            if r["mode"] == "amr_mpi_high"
        }
        for r in sorted(amr_speed, key=lambda rr: (int(rr["N_or_effectiveN"]), int(rr["cores"]))):
            n = int(r["N_or_effectiveN"])
            c = int(r["cores"])
            vals = grouped_amr_mpi.get((str(n), str(c)), [])
            rep = len(vals)
            med = statistics.median(vals) if vals else math.nan
            sd = std(vals)
            lines.append(
                f"| {n} | {c} | {rep} | {f64(r['walltime_sec']):.6f} | {med:.6f} | {sd:.6f} | {f64(r['speedup']):.3f} | {f64(amr_eff[(r['N_or_effectiveN'], c)]):.3f} |"
            )
        lines.append("")

    # Best speedup for highest effective resolution.
    uni_3072 = [f64(r["walltime_sec"]) for r in mpi_rows if i64(r["N"]) == 3072 and i64(r["cores"]) == 1]
    amr_raw_3072 = [f64(r["walltime_sec"]) for r in amr_mpi_rows if i64(r["effective_N"]) == 3072]
    amr_grouped_3072 = {
        int(core): vals
        for (eff_n, core), vals in grouped_amr_mpi.items()
        if int(eff_n) == 3072
    }
    amr_best = min(amr_raw_3072, default=math.nan)
    best_repeated = math.nan
    best_repeated_core = None
    if amr_grouped_3072:
        best_repeated_core, best_repeated_vals = min(
            amr_grouped_3072.items(),
            key=lambda item: mean(item[1]),
        )
        best_repeated = mean(best_repeated_vals)
        best_repeated_std = std(best_repeated_vals)
        best_repeated_n = len(best_repeated_vals)
    if uni_3072 and amr_best == amr_best and amr_best > 0:
        lines.append("## Best High-Resolution AMR vs Uniform")
        lines.append("")
        u1 = mean(uni_3072)
        lines.append(f"- Uniform N=3072, p=1 mean runtime: **{u1:.6f} s**")
        lines.append(f"- Fastest measured AMR runtime (effective N=3072): **{amr_best:.6f} s**")
        lines.append(f"- Fastest measured speed-up vs uniform N=3072 serial: **{u1 / amr_best:.3f}x**")
        if best_repeated == best_repeated and best_repeated > 0:
            lines.append(
                f"- Fastest repeated AMR configuration: **p={best_repeated_core}**, "
                f"**{best_repeated:.6f} \\pm {best_repeated_std:.6f} s** ({best_repeated_n} runs)"
            )
            lines.append(
                f"- Fastest repeated-configuration speed-up vs uniform N=3072 serial: **{u1 / best_repeated:.3f}x**"
            )
        lines.append("")

    summary_path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Summarize Task 3 quadrant timing data and generate plots.")
    ap.add_argument("--results-dir", required=True, help="results_task3_quadrant_YYYYMMDD_HHMM directory")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    ensure_dirs(results_dir)

    time_uniform = read_csv(results_dir / "time_uniform.csv")
    time_mpi = read_csv(results_dir / "time_mpi.csv")
    time_amr_match = read_csv(results_dir / "time_amr_match.csv")
    time_amr_mpi = read_csv(results_dir / "time_amr_mpi.csv")
    sensitivity = read_csv(results_dir / "sensitivity.csv")

    # Accept either legacy (`cores`, `run_id`) or simplified (`p`, `repeat`) MPI schema.
    for r in time_mpi:
        if "cores" not in r and "p" in r:
            r["cores"] = r["p"]
        if "run_id" not in r and "repeat" in r:
            r["run_id"] = f"r{r['repeat']}"
        if "AMR_on" not in r:
            r["AMR_on"] = "0"

    mpi_grouped = aggregate_rows(time_mpi, ["N", "cores"], "walltime_sec")
    amr_mpi_grouped = aggregate_rows(time_amr_mpi, ["effective_N", "cores"], "walltime_sec")
    speed_uni, eff_uni = compute_speedup_efficiency("uniform_mpi", mpi_grouped)
    speed_amr, eff_amr = compute_speedup_efficiency("amr_mpi_high", amr_mpi_grouped)

    speedup_rows = speed_uni + speed_amr
    efficiency_rows = eff_uni + eff_amr

    write_csv(
        results_dir / "speedup.csv",
        ["mode", "N_or_effectiveN", "cores", "walltime_sec", "speedup"],
        speedup_rows,
    )
    write_csv(
        results_dir / "efficiency.csv",
        ["mode", "N_or_effectiveN", "cores", "walltime_sec", "efficiency"],
        efficiency_rows,
    )

    if time_uniform:
        plot_time_vs_n(results_dir, time_uniform)
    if time_mpi or time_amr_mpi:
        plot_time_vs_cores(results_dir, time_mpi, time_amr_mpi)
    if speedup_rows:
        plot_speedup_vs_cores(results_dir, speedup_rows)
    if time_uniform and time_amr_match:
        plot_amr_vs_uniform(results_dir, time_uniform, time_amr_match)
    if sensitivity:
        plot_sensitivity(results_dir, sensitivity)

    write_summary_md(
        results_dir / "tables" / "task3_summary.md",
        time_uniform,
        time_mpi,
        time_amr_match,
        time_amr_mpi,
        speedup_rows,
        efficiency_rows,
    )

    print(f"Wrote summary tables/plots into {results_dir}")


if __name__ == "__main__":
    main()
