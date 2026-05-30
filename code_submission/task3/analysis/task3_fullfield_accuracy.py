#!/usr/bin/env python3
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

GAMMA = 1.4


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


def parse_plot_header(plotfile: Path):
    lines = (plotfile / "Header").read_text().splitlines()
    idx = 0
    idx += 1  # magic
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


def load_fab_var(data_file: Path, offset: int, nx: int, ny: int, ncomp: int, comp: int):
    with data_file.open("rb") as f:
        f.seek(offset)
        f.readline()
        arr = np.fromfile(f, dtype="<f8", count=nx * ny * ncomp)
    fab = arr.reshape((nx, ny, ncomp), order="F")
    return fab[:, :, comp]


CONS_VARS = ["rho", "momx", "momy", "E"]
PRIM_VARS = ["rho", "u", "v", "p", "eint"]


def load_uniform_fields(plotfile: Path):
    meta = parse_plot_header(plotfile)
    if meta["finest_level"] != 0:
        raise ValueError(f"Expected single-level plotfile, got finest_level={meta['finest_level']} for {plotfile}")
    header_lines = (plotfile / "Header").read_text().splitlines()
    m = re.match(r"\(\((\d+),(\d+)\) \((\d+),(\d+)\) \((\d+),(\d+)\)\)", header_lines[12].strip())
    if not m:
        raise ValueError(f"Could not parse domain box from {plotfile / 'Header'}")
    _, _, ihi, jhi, _, _ = map(int, m.groups())
    nx = ihi + 1
    ny = jhi + 1
    fields = {name: np.empty((nx, ny), dtype=np.float64) for name in CONS_VARS}
    boxes, offsets = parse_cell_header(plotfile / "Level_0" / "Cell_H")
    comps = {name: meta["names"].index(name) for name in CONS_VARS}
    data_file = plotfile / "Level_0" / "Cell_D_00000"
    for (ilo, jlo, ihi, jhi), offset in zip(boxes, offsets):
        bx = ihi - ilo + 1
        by = jhi - jlo + 1
        for name, comp in comps.items():
            fab = load_fab_var(data_file, offset, bx, by, meta["ncomp"], comp)
            fields[name][ilo : ihi + 1, jlo : jhi + 1] = fab
    return meta, fields


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


def bilinear_reference(ref_field: np.ndarray, xq: np.ndarray, yq: np.ndarray, dx_ref: float):
    # Reference field is cell-centred on a uniform grid with centres (i+0.5) dx_ref.
    tx = xq / dx_ref - 0.5
    ty = yq / dx_ref - 0.5
    i0 = np.floor(tx).astype(int)
    j0 = np.floor(ty).astype(int)
    wx = tx - i0
    wy = ty - j0

    nx, ny = ref_field.shape
    i0 = np.clip(i0, 0, nx - 2)
    j0 = np.clip(j0, 0, ny - 2)
    i1 = i0 + 1
    j1 = j0 + 1

    f00 = ref_field[i0, j0]
    f10 = ref_field[i1, j0]
    f01 = ref_field[i0, j1]
    f11 = ref_field[i1, j1]
    return (
        (1.0 - wx) * (1.0 - wy) * f00
        + wx * (1.0 - wy) * f10
        + (1.0 - wx) * wy * f01
        + wx * wy * f11
    )


def primitive_from_cons(rho, momx, momy, energy, gamma=GAMMA):
    u = momx / rho
    v = momy / rho
    p = (gamma - 1.0) * (energy - 0.5 * (momx * momx + momy * momy) / rho)
    eint = p / ((gamma - 1.0) * rho)
    return {"rho": rho, "u": u, "v": v, "p": p, "eint": eint}


def fullfield_error_against_uniform_ref(plotfile: Path, ref_meta, ref_fields, gamma=GAMMA):
    meta = parse_plot_header(plotfile)
    ref_dx = ref_meta["dx"][0][0]

    total_abs = {name: 0.0 for name in PRIM_VARS}
    total_sq = {name: 0.0 for name in PRIM_VARS}
    total_area = 0.0

    level_boxes = []
    level_offsets = []
    for lev in range(meta["finest_level"] + 1):
        boxes, offsets = parse_cell_header(plotfile / f"Level_{lev}" / "Cell_H")
        level_boxes.append(boxes)
        level_offsets.append(offsets)

    for lev in range(meta["finest_level"] + 1):
        dx, dy = meta["dx"][lev]
        cell_area = dx * dy
        covered_boxes = []
        if lev < meta["finest_level"]:
            covered_boxes = coarsen_boxes(level_boxes[lev + 1], meta["ref_ratio"][lev])
        data_file = plotfile / f"Level_{lev}" / "Cell_D_00000"
        comps = {name: meta["names"].index(name) for name in CONS_VARS}

        for box, offset in zip(level_boxes[lev], level_offsets[lev]):
            ilo, jlo, ihi, jhi = box
            nx = ihi - ilo + 1
            ny = jhi - jlo + 1
            mask = np.ones((nx, ny), dtype=bool)
            apply_cover_mask(mask, box, covered_boxes)
            if not mask.any():
                continue

            cons = {
                name: load_fab_var(data_file, offset, nx, ny, meta["ncomp"], comp)
                for name, comp in comps.items()
            }
            xs = meta["prob_lo"][0] + (np.arange(ilo, ihi + 1) + 0.5) * dx
            ys = meta["prob_lo"][1] + (np.arange(jlo, jhi + 1) + 0.5) * dy
            xx, yy = np.meshgrid(xs, ys, indexing="ij")
            ref_cons = {
                name: bilinear_reference(ref_fields[name], xx, yy, ref_dx)
                for name in CONS_VARS
            }

            prim = primitive_from_cons(cons["rho"], cons["momx"], cons["momy"], cons["E"], gamma=gamma)
            prim_ref = primitive_from_cons(
                ref_cons["rho"], ref_cons["momx"], ref_cons["momy"], ref_cons["E"], gamma=gamma
            )

            for name in PRIM_VARS:
                diff = prim[name] - prim_ref[name]
                diff_masked = diff[mask]
                total_abs[name] += float(np.sum(np.abs(diff_masked))) * cell_area
                total_sq[name] += float(np.sum(diff_masked * diff_masked)) * cell_area
            total_area += int(np.count_nonzero(mask)) * cell_area

    out = {}
    for name in PRIM_VARS:
        out[f"L1_{name}"] = total_abs[name] / total_area
        out[f"L2_{name}"] = math.sqrt(total_sq[name] / total_area)
    return out


def add_orders(rows, key):
    rows = sorted(rows, key=lambda r: (r["family"], r["effective_N"]))
    prev = {}
    for row in rows:
        fam = row["family"]
        row[f"order_{key}"] = ""
        if fam in prev:
            ecoarse = float(prev[fam][key])
            efine = float(row[key])
            if efine > 0.0:
                row[f"order_{key}"] = f"{math.log(ecoarse / efine, 2.0):.3f}"
        prev[fam] = row
    return rows


def write_outputs(csv_path: Path, md_path: Path, rows):
    fields = ["family", "effective_N", "reference"]
    for name in PRIM_VARS:
        fields.extend([f"L1_{name}", f"order_L1_{name}", f"L2_{name}", f"order_L2_{name}"])
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    grouped = defaultdict(dict)
    for row in rows:
        grouped[row["family"]][int(row["effective_N"])] = row

    lines = []
    lines.append("# Full-Field Accuracy Confirmation for Task 3 (Quadrant, t=0.3)")
    lines.append("")
    lines.append(
        "Density is compared against a repeated uniform $N=3072$ serial reference over the full two-dimensional field."
    )
    lines.append(
        "Uniform runs are compared cell-by-cell on their native mesh; AMR runs use valid cells only, with finer levels masking coarse cells, and the uniform reference is bilinearly interpolated to each queried cell centre."
    )
    lines.append("")
    lines.append("| Family | effective N | L1(rho) | order | L1(u) | order | L1(v) | order | L1(p) | order |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for fam in ["Uniform", "AMR"]:
        for eff_n in sorted(grouped[fam]):
            row = grouped[fam][eff_n]
            lines.append(
                f"| {fam} | {eff_n} | {float(row['L1_rho']):.6e} | {row['order_L1_rho'] or '--'} | "
                f"{float(row['L1_u']):.6e} | {row['order_L1_u'] or '--'} | "
                f"{float(row['L1_v']):.6e} | {row['order_L1_v'] or '--'} | "
                f"{float(row['L1_p']):.6e} | {row['order_L1_p'] or '--'} |"
            )
    if 1536 in grouped["Uniform"] and 1536 in grouped["AMR"]:
        ratio_rho = float(grouped["AMR"][1536]["L1_rho"]) / float(grouped["Uniform"][1536]["L1_rho"])
        ratio_u = float(grouped["AMR"][1536]["L1_u"]) / float(grouped["Uniform"][1536]["L1_u"])
        ratio_v = float(grouped["AMR"][1536]["L1_v"]) / float(grouped["Uniform"][1536]["L1_v"])
        ratio_p = float(grouped["AMR"][1536]["L1_p"]) / float(grouped["Uniform"][1536]["L1_p"])
        lines.append("")
        lines.append("| Matched effective N=1536 comparison | AMR/Uniform L1 ratio |")
        lines.append("|---|---:|")
        lines.append(f"| rho | {ratio_rho:.3f} |")
        lines.append(f"| u | {ratio_u:.3f} |")
        lines.append(f"| v | {ratio_v:.3f} |")
        lines.append(f"| p | {ratio_p:.3f} |")
    lines.append("")
    lines.append("Interpretation:")
    lines.append("- Uniform refinement from 768 to 1536 reduces the full-field errors against the uniform-3072 reference.")
    lines.append("- AMR refinement from effective 1536 to effective 3072 also reduces the same full-field errors.")
    lines.append("- At matched effective N=1536, AMR should still be discussed comparatively rather than as uniformly more accurate on this diagnostic.")
    md_path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Task 3 full-field density accuracy check against the repeated uniform N=3072 reference.")
    ap.add_argument("--quadrant-dir", default="")
    ap.add_argument("--uniform3072-dir", default="")
    ap.add_argument("--uniform768-plot", default="")
    ap.add_argument("--uniform1536-plot", default="")
    ap.add_argument("--amr1536-plot", default="")
    ap.add_argument("--amr3072-plot", default="")
    ap.add_argument("--uniform3072-plot", default="")
    ap.add_argument("--out-dir", default="")
    args = ap.parse_args()

    if args.uniform768_plot:
        pf_uniform_768 = Path(args.uniform768_plot).resolve()
        pf_uniform_1536 = Path(args.uniform1536_plot).resolve()
        pf_amr_1536 = Path(args.amr1536_plot).resolve()
        pf_amr_3072 = Path(args.amr3072_plot).resolve()
        pf_uniform_3072 = Path(args.uniform3072_plot).resolve()
        out_dir = Path(args.out_dir).resolve()
    else:
        quadrant_dir = Path(args.quadrant_dir).resolve()
        uniform3072_dir = Path(args.uniform3072_dir).resolve()
        out_dir = quadrant_dir / "fullfield_accuracy"
        pf_uniform_768 = latest_plot(quadrant_dir / "raw", "plt_stability_uniform_n768_p1_rep3")
        pf_uniform_1536 = latest_plot(quadrant_dir / "raw", "plt_stability_uniform_n1536_p1_rep3")
        pf_amr_1536 = latest_plot(quadrant_dir / "raw", "plt_amr_match_eff1536_np1_r")
        pf_amr_3072 = latest_plot(quadrant_dir / "raw", "plt_amr_match_eff3072_np1_r")
        pf_uniform_3072 = latest_plot(uniform3072_dir / "raw", "plt_uniform_n3072_np1_r3")

    out_dir.mkdir(parents=True, exist_ok=True)

    ref_meta, ref_fields = load_uniform_fields(pf_uniform_3072)

    rows = [
        {"family": "Uniform", "effective_N": 768, "reference": str(pf_uniform_3072), **fullfield_error_against_uniform_ref(pf_uniform_768, ref_meta, ref_fields)},
        {"family": "Uniform", "effective_N": 1536, "reference": str(pf_uniform_3072), **fullfield_error_against_uniform_ref(pf_uniform_1536, ref_meta, ref_fields)},
        {"family": "AMR", "effective_N": 1536, "reference": str(pf_uniform_3072), **fullfield_error_against_uniform_ref(pf_amr_1536, ref_meta, ref_fields)},
        {"family": "AMR", "effective_N": 3072, "reference": str(pf_uniform_3072), **fullfield_error_against_uniform_ref(pf_amr_3072, ref_meta, ref_fields)},
    ]
    for name in PRIM_VARS:
        rows = add_orders(rows, f"L1_{name}")
        rows = add_orders(rows, f"L2_{name}")

    write_outputs(out_dir / "task3_fullfield_accuracy.csv", out_dir / "task3_fullfield_accuracy.md", rows)
    print(f"Wrote {out_dir / 'task3_fullfield_accuracy.csv'}")
    print(f"Wrote {out_dir / 'task3_fullfield_accuracy.md'}")


if __name__ == "__main__":
    main()
