#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import numpy as np


CONS_VARS = ["rho", "momx", "momy", "E"]


def latest_plot(raw_dir: Path, prefix: str) -> Path:
    cands = []
    for p in raw_dir.glob(prefix + "*"):
        suffix = p.name[len(prefix) :]
        if suffix.isdigit():
            cands.append((int(suffix), p))
    if not cands:
        raise FileNotFoundError(f"No plotfiles found for prefix {prefix} in {raw_dir}")
    cands.sort()
    return cands[-1][1]


def parse_plot_header(plotfile: Path):
    lines = (plotfile / "Header").read_text().splitlines()
    idx = 0
    idx += 1
    ncomp = int(lines[idx].strip())
    idx += 1
    names = lines[idx : idx + ncomp]
    idx += ncomp
    idx += 1  # dim
    idx += 1  # time
    finest_level = int(lines[idx].strip())
    idx += 1
    prob_lo = tuple(float(v) for v in lines[idx].split())
    idx += 1
    prob_hi = tuple(float(v) for v in lines[idx].split())
    idx += 1
    ref_ratio = [int(v) for v in lines[idx].split()]
    idx += 1
    idx += 2
    dx = []
    for _ in range(finest_level + 1):
        dx.append(tuple(float(v) for v in lines[idx].split()))
        idx += 1
    return {
        "ncomp": ncomp,
        "names": names,
        "finest_level": finest_level,
        "prob_lo": prob_lo,
        "prob_hi": prob_hi,
        "ref_ratio": ref_ratio,
        "dx": dx,
    }


def parse_cell_header(cell_h_path: Path):
    import re

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


def integrate_plotfile(plotfile: Path):
    meta = parse_plot_header(plotfile)
    level_boxes = []
    level_offsets = []
    for lev in range(meta["finest_level"] + 1):
        boxes, offsets = parse_cell_header(plotfile / f"Level_{lev}" / "Cell_H")
        level_boxes.append(boxes)
        level_offsets.append(offsets)

    comps = {name: meta["names"].index(name) for name in CONS_VARS}
    totals = {name: 0.0 for name in CONS_VARS}

    for lev in range(meta["finest_level"] + 1):
        dx, dy = meta["dx"][lev]
        cell_area = dx * dy
        covered_boxes = []
        if lev < meta["finest_level"]:
            covered_boxes = coarsen_boxes(level_boxes[lev + 1], meta["ref_ratio"][lev])
        data_file = plotfile / f"Level_{lev}" / "Cell_D_00000"
        with data_file.open("rb") as f:
            for box, offset in zip(level_boxes[lev], level_offsets[lev]):
                ilo, jlo, ihi, jhi = box
                nx = ihi - ilo + 1
                ny = jhi - jlo + 1
                mask = np.ones((nx, ny), dtype=bool)
                apply_cover_mask(mask, box, covered_boxes)
                if not mask.any():
                    continue

                f.seek(offset)
                f.readline()
                arr = np.fromfile(f, dtype="<f8", count=nx * ny * meta["ncomp"])
                fab = arr.reshape((nx, ny, meta["ncomp"]), order="F")
                for name, comp in comps.items():
                    vals = fab[:, :, comp][mask]
                    totals[name] += float(np.sum(vals)) * cell_area

    return totals


def rel_err(num, ref):
    if ref == 0.0:
        return 0.0 if num == 0.0 else float("nan")
    return abs(num - ref) / abs(ref)


def write_outputs(out_dir: Path, rows):
    csv_path = out_dir / "smooth_conservation.csv"
    fields = [
        "effective_N",
        "L1_level1_coverage_pct",
        "L2_level2_coverage_pct",
        "rho_ref",
        "rho_num",
        "rho_rel_err",
        "momx_ref",
        "momx_num",
        "momx_rel_err",
        "momy_ref",
        "momy_num",
        "momy_rel_err",
        "E_ref",
        "E_num",
        "E_rel_err",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    lines = []
    lines.append("# Smooth Periodic Conservation Check")
    lines.append("")
    lines.append(
        "For the periodic smooth multilevel AMR runs, final domain-integrated conserved quantities are compared against the exact shifted reference plotfile on the same effective finest mesh."
    )
    lines.append("")
    lines.append("| Effective N | L1 coverage (%) | L2 coverage (%) | rel.err mass | rel.err momx | rel.err momy | rel.err energy |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|")
    for row in rows:
        lines.append(
            f"| {row['effective_N']} | {float(row['L1_level1_coverage_pct']):.2f} | "
            f"{float(row['L2_level2_coverage_pct']):.2f} | "
            f"{float(row['rho_rel_err']):.3e} | {float(row['momx_rel_err']):.3e} | "
            f"{float(row['momy_rel_err']):.3e} | {float(row['E_rel_err']):.3e} |"
        )
    (out_dir / "smooth_conservation.md").write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Conservation check for periodic smooth AMR runs.")
    ap.add_argument("--results-dir", required=True)
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    raw_dir = results_dir / "raw"
    checks_dir = results_dir / "checks"
    checks_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for eff_n in [128, 256, 512]:
        num_plot = latest_plot(raw_dir, f"plt_smooth_amr_eff{eff_n}")
        ref_plot = latest_plot(raw_dir, f"plt_smooth_ref_eff{eff_n}")
        num = integrate_plotfile(num_plot)
        ref = integrate_plotfile(ref_plot)

        # Coverage is available in the existing convergence csv.
        conv_csv = checks_dir / "smooth_amr_convergence.csv"
        l1_cov = 0.0
        l2_cov = 0.0
        if conv_csv.exists():
            with conv_csv.open() as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if int(row["effective_N"]) == eff_n:
                        l1_cov = float(row.get("level1_coverage_pct", 0.0))
                        l2_cov = float(row.get("level2_coverage_pct", 0.0))
                        break

        rows.append(
            {
                "effective_N": eff_n,
                "L1_level1_coverage_pct": l1_cov,
                "L2_level2_coverage_pct": l2_cov,
                "rho_ref": ref["rho"],
                "rho_num": num["rho"],
                "rho_rel_err": rel_err(num["rho"], ref["rho"]),
                "momx_ref": ref["momx"],
                "momx_num": num["momx"],
                "momx_rel_err": rel_err(num["momx"], ref["momx"]),
                "momy_ref": ref["momy"],
                "momy_num": num["momy"],
                "momy_rel_err": rel_err(num["momy"], ref["momy"]),
                "E_ref": ref["E"],
                "E_num": num["E"],
                "E_rel_err": rel_err(num["E"], ref["E"]),
            }
        )

    write_outputs(checks_dir, rows)
    print(f"Wrote {checks_dir / 'smooth_conservation.csv'}")
    print(f"Wrote {checks_dir / 'smooth_conservation.md'}")


if __name__ == "__main__":
    main()
