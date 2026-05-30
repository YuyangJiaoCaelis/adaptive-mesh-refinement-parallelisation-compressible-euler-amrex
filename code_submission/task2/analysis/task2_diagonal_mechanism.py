#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from analyze_task2_2d_fullmatrix import parse_boxes, find_tool, parse_run_index


def refined_fraction(level_boxes, base_n, ref_ratio=2):
    cells = 0
    for ilo, jlo, ihi, jhi in level_boxes:
        cells += (ihi - ilo + 1) * (jhi - jlo + 1)
    fine_n = base_n * ref_ratio
    return cells / float(fine_n * fine_n)


def read_diag_metric(metrics_csv: Path, test: int, mode: str):
    with metrics_csv.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if (
                int(row["test"]) == test
                and row["mode"] == mode
                and row["metric"] == "diag_xy"
                and row["case"] == "diag_xslice_vs_yslice"
                and row["variable"] == "rho"
            ):
                return float(row["L1"]), float(row["L2"]), float(row["Linf"])
    raise RuntimeError(f"Could not find diagonal metric for test={test}, mode={mode}")


def main():
    ap = argparse.ArgumentParser(description="Summarize representative diagonal AMR fragmentation for Task 2.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--test", type=int, default=4)
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    checks_dir = results_dir / "checks"
    checks_dir.mkdir(parents=True, exist_ok=True)

    run_index = parse_run_index(results_dir / "logs" / "run_index.tsv")
    run_map = {(r["test"], r["ic"], r["mode"]): r for r in run_index}
    fboxinfo = find_tool("fboxinfo")

    rows = []
    for label, ic in [("x-split", 0), ("diagonal", 2)]:
        case = run_map[(args.test, ic, "amr")]
        boxes = parse_boxes(case["plotfile"], fboxinfo)
        lev1 = boxes.get(1, [])
        rows.append(
            {
                "test": args.test,
                "case": label,
                "plotfile": str(case["plotfile"]),
                "level1_boxes": len(lev1),
                "level1_refined_fraction": refined_fraction(lev1, base_n=200),
                "level1_refined_percent": 100.0 * refined_fraction(lev1, base_n=200),
            }
        )

    uni_l1, _, _ = read_diag_metric(checks_dir / "task2_metrics_long.csv", args.test, "uniform")
    amr_l1, _, _ = read_diag_metric(checks_dir / "task2_metrics_long.csv", args.test, "amr")
    ratio = amr_l1 / uni_l1 if uni_l1 > 0.0 else float("nan")

    csv_path = checks_dir / "task2_diagonal_mechanism.csv"
    with csv_path.open("w", newline="") as f:
        fieldnames = [
            "test",
            "case",
            "plotfile",
            "level1_boxes",
            "level1_refined_fraction",
            "level1_refined_percent",
        ]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    lines = []
    lines.append(f"# Task 2 Diagonal Mechanism Check (Test {args.test})")
    lines.append("")
    lines.append("| Case | Level-1 box count | Level-1 refined area |")
    lines.append("|---|---:|---:|")
    for row in rows:
        lines.append(
            f"| {row['case']} | {row['level1_boxes']} | {row['level1_refined_percent']:.2f}\\% |"
        )
    lines.append("")
    lines.append(
        f"Diagonal symmetry metric in `checks/task2_metrics_long.csv`: uniform L1 = {uni_l1:.6e}, "
        f"AMR L1 = {amr_l1:.6e}, ratio = {ratio:.3f}."
    )
    (checks_dir / "task2_diagonal_mechanism.md").write_text("\n".join(lines) + "\n")

    print(f"CSV={csv_path}")
    print(f"MD={checks_dir / 'task2_diagonal_mechanism.md'}")


if __name__ == "__main__":
    main()
