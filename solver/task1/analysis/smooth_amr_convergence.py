#!/usr/bin/env python3
import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


VARS = ["rho", "momx", "momy", "E"]


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


def parse_pressure_floor(log_path: Path):
    text = log_path.read_text()
    match = re.search(
        r"Diag pressure-floor triggers \(sum over MPI ranks\): stage1=(\d+) stage2=(\d+) total=(\d+)",
        text,
    )
    if match:
        return tuple(int(v) for v in match.groups())
    return (0, 0, 0)


def parse_runtime(log_path: Path):
    text = log_path.read_text()
    match = re.search(r"Run time =\s*([0-9.eE+-]+)", text)
    return float(match.group(1)) if match else math.nan


def parse_plot_header(plotfile: Path):
    lines = (plotfile / "Header").read_text().splitlines()
    idx = 0
    _magic = lines[idx]
    idx += 1
    ncomp = int(lines[idx].strip())
    idx += 1
    names = lines[idx : idx + ncomp]
    idx += ncomp
    dim = int(lines[idx].strip())
    idx += 1
    _time = float(lines[idx].strip())
    idx += 1
    finest_level = int(lines[idx].strip())
    idx += 1
    prob_lo = tuple(float(v) for v in lines[idx].split())
    idx += 1
    prob_hi = tuple(float(v) for v in lines[idx].split())
    idx += 1
    ref_ratio = [int(v) for v in lines[idx].split()]
    idx += 1
    idx += 2  # domain boxes and level steps
    dx = []
    for _ in range(finest_level + 1):
        dx.append(tuple(float(v) for v in lines[idx].split()))
        idx += 1
    return {
        "ncomp": ncomp,
        "names": names,
        "dim": dim,
        "finest_level": finest_level,
        "prob_lo": prob_lo,
        "prob_hi": prob_hi,
        "ref_ratio": ref_ratio,
        "dx": dx,
    }


def parse_cell_header(cell_h_path: Path):
    text = cell_h_path.read_text().splitlines()
    boxes = []
    offsets = []
    box_pat = re.compile(r"\(\((\d+),(\d+)\) \((\d+),(\d+)\) \((\d+),(\d+)\)\)")
    for line in text:
        s = line.strip()
        m = box_pat.match(s)
        if m:
            ilo, jlo, ihi, jhi, _, _ = map(int, m.groups())
            boxes.append((ilo, jlo, ihi, jhi))
        elif s.startswith("FabOnDisk:"):
            offsets.append(int(s.split()[-1]))
    if len(boxes) != len(offsets):
        raise ValueError(f"Mismatched box/offset count in {cell_h_path}")
    return boxes, offsets


def coarsen_boxes(boxes, ref_ratio):
    out = []
    for ilo, jlo, ihi, jhi in boxes:
        out.append((ilo // ref_ratio, jlo // ref_ratio, ihi // ref_ratio, jhi // ref_ratio))
    return out


def apply_cover_mask(mask, box, covered_boxes):
    ilo, jlo, ihi, jhi = box
    for cilo, cjlo, cihi, cjhi in covered_boxes:
        oi0 = max(ilo, cilo)
        oj0 = max(jlo, cjlo)
        oi1 = min(ihi, cihi)
        oj1 = min(jhi, cjhi)
        if oi0 <= oi1 and oj0 <= oj1:
            mask[(oi0 - ilo) : (oi1 - ilo + 1), (oj0 - jlo) : (oj1 - jlo + 1)] = False


def exact_state(x, y, *, x0, y0, rho0, amp, u0, v0, p0, sigma, profile, kx, ky):
    if profile == "gaussian":
        r2 = (x - x0) ** 2 + (y - y0) ** 2
        rho = rho0 + amp * np.exp(-r2 / (sigma * sigma))
    elif profile == "sinusoid":
        rho = rho0 + amp * np.sin(2.0 * math.pi * (kx * (x - x0) + ky * (y - y0)))
    else:
        raise ValueError(f"Unknown smooth profile: {profile}")
    momx = rho * u0
    momy = rho * v0
    energy = p0 / (1.4 - 1.0) + 0.5 * rho * (u0 * u0 + v0 * v0)
    return {
        "rho": rho,
        "momx": momx,
        "momy": momy,
        "E": energy,
    }


def compute_amr_errors(plotfile: Path, params):
    meta = parse_plot_header(plotfile)
    level_boxes = []
    level_offsets = []
    for lev in range(meta["finest_level"] + 1):
        boxes, offsets = parse_cell_header(plotfile / f"Level_{lev}" / "Cell_H")
        level_boxes.append(boxes)
        level_offsets.append(offsets)

    accum_abs = {v: 0.0 for v in VARS}
    accum_sq = {v: 0.0 for v in VARS}
    volume = 0.0
    coverage_pct = []
    domain_area = (meta["prob_hi"][0] - meta["prob_lo"][0]) * (meta["prob_hi"][1] - meta["prob_lo"][1])

    for lev in range(meta["finest_level"] + 1):
        dx, dy = meta["dx"][lev]
        cell_area = dx * dy
        covered_boxes = []
        if lev < meta["finest_level"]:
            covered_boxes = coarsen_boxes(level_boxes[lev + 1], meta["ref_ratio"][lev])

        level_valid_area = 0.0
        data_file = plotfile / f"Level_{lev}" / "Cell_D_00000"
        with data_file.open("rb") as f:
            for box, offset in zip(level_boxes[lev], level_offsets[lev]):
                ilo, jlo, ihi, jhi = box
                nx = ihi - ilo + 1
                ny = jhi - jlo + 1
                f.seek(offset)
                f.readline()
                arr = np.fromfile(f, dtype="<f8", count=nx * ny * meta["ncomp"])
                fab = arr.reshape((nx, ny, meta["ncomp"]), order="F")

                mask = np.ones((nx, ny), dtype=bool)
                apply_cover_mask(mask, box, covered_boxes)
                if not mask.any():
                    continue

                xs = meta["prob_lo"][0] + (np.arange(ilo, ihi + 1) + 0.5) * dx
                ys = meta["prob_lo"][1] + (np.arange(jlo, jhi + 1) + 0.5) * dy
                xx, yy = np.meshgrid(xs, ys, indexing="ij")
                exact = exact_state(xx, yy, **params)

                for n, name in enumerate(meta["names"]):
                    diff = fab[:, :, n] - exact[name]
                    diff_masked = diff[mask]
                    accum_abs[name] += float(np.sum(np.abs(diff_masked))) * cell_area
                    accum_sq[name] += float(np.sum(diff_masked * diff_masked)) * cell_area

                valid_cells = int(np.count_nonzero(mask))
                level_valid_area += valid_cells * cell_area
                volume += valid_cells * cell_area

        coverage_pct.append(100.0 * level_valid_area / domain_area)

    rows = {}
    for name in VARS:
        rows[f"L1_{name}"] = accum_abs[name] / volume
        rows[f"L2_{name}"] = math.sqrt(accum_sq[name] / volume)
    return rows, coverage_pct


def add_orders(rows, key):
    rows = sorted(rows, key=lambda r: r["effective_N"])
    prev = None
    for row in rows:
        row[f"order_{key}"] = ""
        if prev is not None:
            ecoarse = float(prev[key])
            efine = float(row[key])
            if efine > 0.0:
                row[f"order_{key}"] = f"{math.log(ecoarse / efine, 2.0):.3f}"
        prev = row
    return rows


def write_markdown(md_path: Path, rows, *, profile, max_level):
    lines = []
    lines.append("# Smooth AMR Convergence")
    lines.append("")
    lines.append(
        f"A {max_level}-level AMR smooth-{profile} advection test is compared against the exact shifted state on the same effective finest mesh."
    )
    lines.append("")
    lines.append("| Effective N | Base N | Level-1 coverage (%) | Level-2 coverage (%) | L1(rho) | order | L2(rho) | order | Pressure-floor triggers | runtime (s) |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| {row['effective_N']} | {row['base_N']} | {float(row['level1_coverage_pct']):.2f} | "
            f"{float(row['level2_coverage_pct']):.2f} | "
            f"{float(row['L1_rho']):.6e} | {row['order_L1_rho'] or '--'} | "
            f"{float(row['L2_rho']):.6e} | {row['order_L2_rho'] or '--'} | "
            f"{row['pressure_floor_total']} | {float(row['runtime_s']):.3f} |"
        )
    md_path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Summarize smooth AMR convergence runs.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--stop-time", type=float, required=True)
    ap.add_argument("--x0", type=float, required=True)
    ap.add_argument("--y0", type=float, required=True)
    ap.add_argument("--rho0", type=float, required=True)
    ap.add_argument("--amp", type=float, required=True)
    ap.add_argument("--u0", type=float, required=True)
    ap.add_argument("--v0", type=float, required=True)
    ap.add_argument("--p0", type=float, required=True)
    ap.add_argument("--sigma", type=float, required=True)
    ap.add_argument("--profile", choices=["gaussian", "sinusoid"], default="gaussian")
    ap.add_argument("--kx", type=float, default=1.0)
    ap.add_argument("--ky", type=float, default=1.0)
    ap.add_argument("--max-level", type=int, default=1)
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    raw_dir = results_dir / "raw"
    log_dir = results_dir / "logs"
    checks_dir = results_dir / "checks"
    checks_dir.mkdir(parents=True, exist_ok=True)

    params = {
        "x0": args.x0 + args.u0 * args.stop_time,
        "y0": args.y0 + args.v0 * args.stop_time,
        "rho0": args.rho0,
        "amp": args.amp,
        "u0": args.u0,
        "v0": args.v0,
        "p0": args.p0,
        "sigma": args.sigma,
        "profile": args.profile,
        "kx": args.kx,
        "ky": args.ky,
    }

    rows = []
    for eff_n in [128, 256, 512]:
        plotfile = latest_plot(raw_dir, f"plt_smooth_amr_eff{eff_n}")
        errs, coverage = compute_amr_errors(plotfile, params)
        stage1, stage2, total = parse_pressure_floor(log_dir / f"smooth_amr_eff{eff_n}.log")
        rows.append(
            {
                "effective_N": eff_n,
                "base_N": eff_n // (2 ** args.max_level),
                "level0_coverage_pct": coverage[0],
                "level1_coverage_pct": coverage[1] if len(coverage) > 1 else 0.0,
                "level2_coverage_pct": coverage[2] if len(coverage) > 2 else 0.0,
                "runtime_s": parse_runtime(log_dir / f"smooth_amr_eff{eff_n}.log"),
                "pressure_floor_stage1": stage1,
                "pressure_floor_stage2": stage2,
                "pressure_floor_total": total,
                **errs,
            }
        )

    rows = add_orders(rows, "L1_rho")
    rows = add_orders(rows, "L2_rho")

    csv_path = checks_dir / "smooth_amr_convergence.csv"
    fieldnames = [
        "effective_N",
        "base_N",
        "level0_coverage_pct",
        "level1_coverage_pct",
        "level2_coverage_pct",
        "runtime_s",
        "pressure_floor_stage1",
        "pressure_floor_stage2",
        "pressure_floor_total",
        "L1_rho",
        "order_L1_rho",
        "L2_rho",
        "order_L2_rho",
        "L1_momx",
        "L2_momx",
        "L1_momy",
        "L2_momy",
        "L1_E",
        "L2_E",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    write_markdown(checks_dir / "smooth_amr_convergence.md", rows, profile=args.profile, max_level=args.max_level)
    print(f"CSV={csv_path}")
    print(f"MD={checks_dir / 'smooth_amr_convergence.md'}")


if __name__ == "__main__":
    main()
