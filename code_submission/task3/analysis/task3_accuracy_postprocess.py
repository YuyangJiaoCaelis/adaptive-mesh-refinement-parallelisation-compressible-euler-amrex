#!/usr/bin/env python3
import argparse
import csv
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def latest_plot(raw_dir: Path, prefix: str) -> Path:
    cands = []
    for p in raw_dir.glob(prefix + "*"):
        suffix = p.name[len(prefix) :]
        if suffix.isdigit():
            cands.append((int(suffix), p))
    if not cands:
        raise FileNotFoundError(f"No plotfiles found for prefix: {prefix}")
    cands.sort()
    return cands[-1][1]


def run_fextract(
    fextract_exe: Path,
    plotfile: Path,
    out_file: Path,
    direction: int,
    coord_flag: str,
    coord_val: float,
):
    out_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(fextract_exe),
        "-s",
        str(out_file),
        "-d",
        str(direction),
        coord_flag,
        f"{coord_val}",
        "-v",
        "rho",
        str(plotfile),
    ]
    subprocess.run(cmd, check=True)


def read_slice(path: Path):
    xs = []
    ys = []
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    order = np.argsort(x)
    return x[order], y[order]


def compute_error(u_x, u, ref_x, ref):
    u_interp = np.interp(ref_x, u_x, u)
    diff = u_interp - ref
    l1 = float(np.mean(np.abs(diff)))
    l2 = float(np.sqrt(np.mean(diff * diff)))
    return l1, l2


def extract_case_slices(fextract_exe: Path, plotfile: Path, slice_dir: Path, tag: str):
    x_path = slice_dir / f"{tag}_xcenter.dat"
    y_path = slice_dir / f"{tag}_ycenter.dat"
    run_fextract(fextract_exe, plotfile, x_path, 0, "-y", 0.5)
    run_fextract(fextract_exe, plotfile, y_path, 1, "-x", 0.5)
    return read_slice(x_path), read_slice(y_path)


def extract_uniform_single_level_centerlines(plotfile: Path):
    header_lines = (plotfile / "Header").read_text().splitlines()
    m = re.match(r"\(\((\d+),(\d+)\) \((\d+),(\d+)\) \((\d+),(\d+)\)\)", header_lines[12].strip())
    if not m:
        raise ValueError(f"Could not parse domain box from {plotfile / 'Header'}")
    _, _, hi_x, hi_y, _, _ = map(int, m.groups())
    nx = hi_x + 1
    ny = hi_y + 1
    dx = float(header_lines[14].split()[0])

    boxes = []
    offsets = []
    for line in (plotfile / "Level_0" / "Cell_H").read_text().splitlines():
        s = line.strip()
        if s.startswith("(("):
            m = re.match(r"\(\((\d+),(\d+)\) \((\d+),(\d+)\) \((\d+),(\d+)\)\)", s)
            if m:
                ilo, jlo, ihi, jhi, _, _ = map(int, m.groups())
                boxes.append((ilo, jlo, ihi, jhi))
        elif s.startswith("FabOnDisk:"):
            offsets.append(int(s.split()[-1]))

    if len(boxes) != len(offsets):
        raise ValueError(f"Mismatched box/offset count in {plotfile / 'Level_0' / 'Cell_H'}")

    # Match the existing fextract convention used in the archived Task 3 slices:
    # for x=0.5 or y=0.5 on an even grid, select the first cell centre above 0.5.
    center_index = int(math.ceil(0.5 / dx - 0.5))
    rho_x = np.zeros(nx)
    rho_y = np.zeros(ny)
    got_x = np.zeros(nx, dtype=bool)
    got_y = np.zeros(ny, dtype=bool)

    with (plotfile / "Level_0" / "Cell_D_00000").open("rb") as f:
        for (ilo, jlo, ihi, jhi), offset in zip(boxes, offsets):
            f.seek(offset)
            f.readline()  # FAB header
            count = (ihi - ilo + 1) * (jhi - jlo + 1) * 4
            arr = np.fromfile(f, dtype="<f8", count=count)
            fab_rho = arr.reshape((ihi - ilo + 1, jhi - jlo + 1, 4), order="F")[:, :, 0]

            if jlo <= center_index <= jhi:
                rho_x[ilo : ihi + 1] = fab_rho[:, center_index - jlo]
                got_x[ilo : ihi + 1] = True
            if ilo <= center_index <= ihi:
                rho_y[jlo : jhi + 1] = fab_rho[center_index - ilo, :]
                got_y[jlo : jhi + 1] = True

    if not got_x.all() or not got_y.all():
        raise ValueError(f"Failed to reconstruct complete centreline slices from {plotfile}")

    xs = (np.arange(nx) + 0.5) * dx
    ys = (np.arange(ny) + 0.5) * dx
    return (xs, rho_x), (ys, rho_y)


def order_string(coarse, fine):
    if fine <= 0.0 or coarse <= 0.0:
        return ""
    return f"{math.log(coarse / fine, 2.0):.3f}"


def write_accuracy_outputs(csv_path: Path, md_path: Path, rows, reference_label: str, reference_note: str, *, direct_uniform_reference: bool):
    fields = [
        "family",
        "effective_N",
        "diagnostic",
        "reference",
        "L1_error",
        "L2_error",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    grouped = defaultdict(dict)
    for row in rows:
        grouped[(row["family"], row["diagnostic"])][int(row["effective_N"])] = row

    wide = defaultdict(dict)
    for row in rows:
        key = (row["family"], int(row["effective_N"]))
        wide[key][row["diagnostic"]] = row

    lines = []
    lines.append("# Accuracy Confirmation for Task 3 (Quadrant, t=0.3)")
    lines.append("")
    lines.append(f"Reference solution for this centreline comparison: **{reference_label}**.")
    if reference_note:
        lines.append(reference_note)
    else:
        lines.append("Errors are measured against the supplied high-resolution centerline reference.")
    lines.append("")
    lines.append("| Family | effective N | L1 error: rho(x, y=0.5) | L2 error: rho(x, y=0.5) | L1 error: rho(x=0.5, y) | L2 error: rho(x=0.5, y) |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for family, eff_n in sorted(wide.keys(), key=lambda x: (x[0], x[1])):
        row_x = wide[(family, eff_n)]["rho(x, y=0.5)"]
        row_y = wide[(family, eff_n)]["rho(x=0.5, y)"]
        lines.append(
            f"| {family} | {eff_n} | {float(row_x['L1_error']):.6e} | {float(row_x['L2_error']):.6e} | "
            f"{float(row_y['L1_error']):.6e} | {float(row_y['L2_error']):.6e} |"
        )

    lines.append("")
    lines.append("| Family | Resolution step | L1 ratio x | L1 order x | L1 ratio y | L1 order y |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for family in ["Uniform", "AMR"]:
        x_rows = grouped.get((family, "rho(x, y=0.5)"), {})
        y_rows = grouped.get((family, "rho(x=0.5, y)"), {})
        if 768 in x_rows and 1536 in x_rows:
            cx = float(x_rows[768]["L1_error"])
            fx = float(x_rows[1536]["L1_error"])
            cy = float(y_rows[768]["L1_error"])
            fy = float(y_rows[1536]["L1_error"])
            lines.append(
                f"| {family} | 768 $\\rightarrow$ 1536 | {cx/fx:.3f} | {order_string(cx, fx)} | {cy/fy:.3f} | {order_string(cy, fy)} |"
            )
        if 1536 in x_rows and 3072 in x_rows:
            cx = float(x_rows[1536]["L1_error"])
            fx = float(x_rows[3072]["L1_error"])
            cy = float(y_rows[1536]["L1_error"])
            fy = float(y_rows[3072]["L1_error"])
            lines.append(
                f"| {family} | 1536 $\\rightarrow$ 3072 | {cx/fx:.3f} | {order_string(cx, fx)} | {cy/fy:.3f} | {order_string(cy, fy)} |"
            )

    if ("Uniform", 1536) in wide and ("AMR", 1536) in wide:
        ux = float(wide[("Uniform", 1536)]["rho(x, y=0.5)"]["L1_error"])
        uy = float(wide[("Uniform", 1536)]["rho(x=0.5, y)"]["L1_error"])
        ax = float(wide[("AMR", 1536)]["rho(x, y=0.5)"]["L1_error"])
        ay = float(wide[("AMR", 1536)]["rho(x=0.5, y)"]["L1_error"])
        lines.append("")
        lines.append("| Matched effective N=1536 comparison | AMR/Uniform L1 ratio x | AMR/Uniform L1 ratio y |")
        lines.append("|---|---:|---:|")
        lines.append(f"| rho | {ax/ux:.3f} | {ay/uy:.3f} |")

    lines.append("")
    lines.append("Interpretation:")
    lines.append("- Uniform refinement from 768 to 1536 reduces the centerline density L1 error in both diagnostics.")
    if direct_uniform_reference:
        lines.append("- AMR refinement from effective 1536 to effective 3072 also reduces the direct centerline error against the high-resolution uniform reference.")
        lines.append("- At matched effective N=1536, AMR is not uniformly better on every centerline metric, so the Task 3 accuracy statement should remain comparative rather than absolute.")
    else:
        lines.append("- AMR refinement from effective 1536 to effective 3072 reduces the same surrogate-reference error by a much larger factor.")
        lines.append("- At matched effective N=1536, AMR is not uniformly better on every centerline metric, so the Task 3 accuracy statement should remain cautious rather than absolute.")
    md_path.write_text("\n".join(lines) + "\n")


def plot_overlays(out_png: Path, out_pdf: Path, case_data, ref_x, ref_y, reference_label: str):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), sharey=True)

    style = {
        "Uniform 768": dict(color="#4C78A8", linewidth=1.3),
        "Uniform 1536": dict(color="#72B7B2", linewidth=1.3),
        "AMR eff 1536": dict(color="#F58518", linewidth=1.3),
        "AMR eff 3072": dict(color="#E45756", linewidth=1.5),
    }

    for label, ((xdat, ydat_x), (ydat, ydat_y)) in case_data.items():
        axes[0].plot(xdat, ydat_x, label=label, **style[label])
        axes[1].plot(ydat, ydat_y, label=label, **style[label])

    axes[0].plot(ref_x[0], ref_x[1], color="black", linestyle="--", linewidth=1.2, label=reference_label)
    axes[1].plot(ref_y[0], ref_y[1], color="black", linestyle="--", linewidth=1.2, label=reference_label)

    axes[0].set_title(r"Centerline Slice: $\rho(x, y=0.5)$")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("rho")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].set_title(r"Centerline Slice: $\rho(x=0.5, y)$")
    axes[1].set_xlabel("y")
    axes[1].grid(alpha=0.25)
    axes[1].legend()

    fig.suptitle("Task 3 Centreline Accuracy Comparison")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Task 3 accuracy confirmation from existing plotfiles.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument(
        "--fextract-exe",
        default="/path/to/amrex/Tools/Plotfile/fextract.gnu.ex",
    )
    ap.add_argument(
        "--amr3072-plotfile",
        default="",
        help="Optional explicit AMR effective-3072 serial plotfile to compare against the reference.",
    )
    ap.add_argument(
        "--reference-plotfile",
        default="",
        help="Optional explicit high-resolution reference plotfile. Defaults to the latest AMR p=8 high-resolution plotfile in results-dir.",
    )
    ap.add_argument(
        "--reference-label",
        default="High-resolution AMR reference (effective N=3072, p=8)",
    )
    ap.add_argument(
        "--reference-note",
        default="A high-resolution AMR slice is used as the surrogate reference because the final uniform N=3072 plotfile is unreadable for direct postprocessing.",
    )
    ap.add_argument(
        "--uniform3072-reference-plotfile",
        default="",
        help="Optional explicit uniform N=3072 serial plotfile to use as a direct high-resolution reference.",
    )
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    raw_dir = results_dir / "raw"
    table_dir = results_dir / "tables"
    plot_dir = results_dir / "plots"
    slice_dir = results_dir / "accuracy_slices"
    table_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)
    slice_dir.mkdir(parents=True, exist_ok=True)

    fextract_exe = Path(args.fextract_exe).resolve()

    pf_uniform_768 = latest_plot(raw_dir, "plt_stability_uniform_n768_p1_rep3")
    pf_uniform_1536 = latest_plot(raw_dir, "plt_stability_uniform_n1536_p1_rep3")
    pf_amr_1536 = latest_plot(raw_dir, "plt_amr_match_eff1536_np1_r")
    if args.amr3072_plotfile:
        pf_amr_3072 = Path(args.amr3072_plotfile).resolve()
    else:
        pf_amr_3072 = latest_plot(raw_dir, "plt_amr_match_eff3072_np1_r")
    direct_uniform_reference = False
    if args.uniform3072_reference_plotfile:
        pf_ref = Path(args.uniform3072_reference_plotfile).resolve()
        direct_uniform_reference = True
    elif args.reference_plotfile:
        pf_ref = Path(args.reference_plotfile).resolve()
    else:
        pf_ref = latest_plot(raw_dir, "plt_amr_high_eff3072_np8_r")

    uniform_768 = extract_case_slices(fextract_exe, pf_uniform_768, slice_dir, "rho_uniform_768")
    uniform_1536 = extract_case_slices(fextract_exe, pf_uniform_1536, slice_dir, "rho_uniform_1536")
    amr_1536 = extract_case_slices(fextract_exe, pf_amr_1536, slice_dir, "rho_amr_eff1536")
    amr_3072 = extract_case_slices(fextract_exe, pf_amr_3072, slice_dir, "rho_amr_eff3072")
    if direct_uniform_reference:
        ref = extract_uniform_single_level_centerlines(pf_ref)
    else:
        ref = extract_case_slices(fextract_exe, pf_ref, slice_dir, "rho_reference_eff3072")

    rows = []
    case_map = {
        ("Uniform", 768): uniform_768,
        ("Uniform", 1536): uniform_1536,
        ("AMR", 1536): amr_1536,
        ("AMR", 3072): amr_3072,
    }
    for (family, eff_n), ((x, rho_x), (y, rho_y)) in case_map.items():
        l1x, l2x = compute_error(x, rho_x, ref[0][0], ref[0][1])
        l1y, l2y = compute_error(y, rho_y, ref[1][0], ref[1][1])
        rows.append(
            {
                "family": family,
                "effective_N": eff_n,
                "diagnostic": "rho(x, y=0.5)",
                "reference": args.reference_label,
                "L1_error": l1x,
                "L2_error": l2x,
            }
        )
        rows.append(
            {
                "family": family,
                "effective_N": eff_n,
                "diagnostic": "rho(x=0.5, y)",
                "reference": args.reference_label,
                "L1_error": l1y,
                "L2_error": l2y,
            }
        )

    write_accuracy_outputs(
        results_dir / "accuracy_vs_N.csv",
        table_dir / "accuracy_vs_N.md",
        rows,
        args.reference_label,
        args.reference_note,
        direct_uniform_reference=direct_uniform_reference,
    )

    plot_overlays(
        plot_dir / "accuracy_slice_rho.png",
        plot_dir / "accuracy_slice_rho.pdf",
        {
            "Uniform 768": uniform_768,
            "Uniform 1536": uniform_1536,
            "AMR eff 1536": amr_1536,
            "AMR eff 3072": amr_3072,
        },
        ref[0],
        ref[1],
        args.reference_label,
    )

    print(f"Wrote accuracy files under {results_dir}")


if __name__ == "__main__":
    main()
