#!/usr/bin/env python3
import argparse
import csv
import math
import re
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


IC_LABEL = {0: "xsplit", 1: "ysplit", 2: "diag"}
IC_TITLE = {0: "x-split", 1: "y-split", 2: "diagonal split"}


def find_tool(name: str) -> Path:
    env_home = Path(__import__("os").environ.get("AMREX_HOME", ""))
    if env_home:
        cand = env_home / "Tools" / "Plotfile" / f"{name}.gnu.ex"
        if cand.exists():
            return cand
    raise FileNotFoundError(f"Could not find {name}.gnu.ex under AMREX_HOME/Tools/Plotfile")


def parse_run_index(path: Path):
    rows = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "test": int(row["test"]),
                    "ic": int(row["ic"]),
                    "split": row["split"],
                    "mode": row["mode"],
                    "plotfile": Path(row["plotfile"]),
                    "log": Path(row["log"]),
                }
            )
    return rows


def parse_plot_geom(plotfile: Path):
    header = plotfile / "Header"
    lines = [ln.strip() for ln in header.read_text().splitlines() if ln.strip()]
    idx = 0
    _version = lines[idx]
    idx += 1
    nvars = int(lines[idx])
    idx += 1
    idx += nvars
    dim = int(lines[idx])
    idx += 1
    time = float(lines[idx])
    idx += 1
    finest = int(lines[idx])
    idx += 1
    prob_lo = [float(x) for x in lines[idx].split()]
    idx += 1
    prob_hi = [float(x) for x in lines[idx].split()]
    idx += 1

    dom_line = None
    for line in lines[idx:]:
        if "((" in line and "))" in line:
            dom_line = line
            break
    if dom_line is None:
        raise RuntimeError(f"Could not parse domain tuple in {header}")

    nums = [int(v) for v in re.findall(r"-?\d+", dom_line)]
    dom_lo = nums[:dim]
    dom_hi = nums[dim : 2 * dim]
    ncell = [dom_hi[d] - dom_lo[d] + 1 for d in range(dim)]
    return {
        "dim": dim,
        "time": time,
        "finest": finest,
        "prob_lo": prob_lo,
        "prob_hi": prob_hi,
        "ncell": ncell,
    }


def parse_fextrema_rho(plotfile: Path, fextrema: Path):
    txt = subprocess.check_output([str(fextrema), "-v", "rho", str(plotfile)], text=True)
    for raw in txt.splitlines():
        line = raw.strip()
        if not line.startswith("rho"):
            continue
        parts = line.split()
        if len(parts) >= 3:
            return float(parts[1]), float(parts[2])
    raise RuntimeError(f"rho extrema not found for {plotfile}")


def parse_boxes(plotfile: Path, fboxinfo: Path):
    txt = subprocess.check_output([str(fboxinfo), "-f", str(plotfile)], text=True)
    levels = {}
    cur_level = None
    lev_re = re.compile(r"^\s*level\s+(\d+)\s*$")
    box_re = re.compile(r"box\s+\d+:\s*\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)\s*\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)")
    for raw in txt.splitlines():
        mlev = lev_re.match(raw)
        if mlev:
            cur_level = int(mlev.group(1))
            levels.setdefault(cur_level, [])
            continue
        mbox = box_re.search(raw)
        if mbox and cur_level is not None:
            levels[cur_level].append(tuple(int(mbox.group(i)) for i in range(1, 5)))
    return levels


def read_ppm_p6(path: Path):
    with path.open("rb") as f:
        def next_token():
            while True:
                ch = f.read(1)
                if not ch:
                    return None
                if ch in b" \t\r\n":
                    continue
                if ch == b"#":
                    f.readline()
                    continue
                token = bytearray(ch)
                break
            while True:
                ch = f.read(1)
                if not ch or ch in b" \t\r\n":
                    break
                if ch == b"#":
                    f.readline()
                    break
                token.extend(ch)
            return bytes(token)

        magic = next_token()
        if magic != b"P6":
            raise RuntimeError(f"Expected P6 PPM, got {magic!r} in {path}")
        nx = int(next_token())
        ny = int(next_token())
        maxv = int(next_token())
        if maxv != 255:
            raise RuntimeError(f"Unsupported maxv={maxv} in {path}")
        raw = f.read(nx * ny * 3)
        if len(raw) != nx * ny * 3:
            raise RuntimeError(f"Short PPM payload in {path}")
    arr = np.frombuffer(raw, dtype=np.uint8).reshape((ny, nx, 3))
    return arr


def make_snapshot(plotfile: Path, vmin: float, vmax: float, out_ppm: Path, fsnapshot: Path):
    subprocess.run(
        [
            str(fsnapshot),
            "-v",
            "rho",
            "-m",
            f"{vmin:.16g}",
            "-M",
            f"{vmax:.16g}",
            str(plotfile),
        ],
        check=True,
    )
    src = Path(f"{plotfile}.rho.ppm")
    if not src.exists():
        src = Path(f"{plotfile.name}.rho.ppm")
    if not src.exists():
        raise FileNotFoundError(f"Could not find fsnapshot output for {plotfile}")
    out_ppm.parent.mkdir(parents=True, exist_ok=True)
    src.replace(out_ppm)


def read_csv_strip(path: Path):
    rows = []
    with path.open() as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        for row in reader:
            if not row:
                continue
            rows.append([x.strip() for x in row])
    return header, rows


def compress_duplicate_x(x: np.ndarray, values: np.ndarray):
    order = np.argsort(x)
    x = x[order]
    values = values[order]
    ux, inv = np.unique(x, return_inverse=True)
    sums = np.zeros_like(ux, dtype=float)
    cnts = np.zeros_like(ux, dtype=float)
    np.add.at(sums, inv, values)
    np.add.at(cnts, inv, 1.0)
    return ux, sums / np.maximum(cnts, 1.0)


def parse_slice(path: Path, gamma: float = 1.4):
    _, rows = read_csv_strip(path)
    x = np.array([float(r[0]) for r in rows], dtype=float)
    rho = np.array([float(r[1]) for r in rows], dtype=float)
    momx = np.array([float(r[2]) for r in rows], dtype=float)
    momy = np.array([float(r[3]) for r in rows], dtype=float)
    ener = np.array([float(r[4]) for r in rows], dtype=float)

    u = momx / rho
    v = momy / rho
    kin = 0.5 * (momx * momx + momy * momy) / rho
    p = (gamma - 1.0) * (ener - kin)
    eint = p / ((gamma - 1.0) * rho)

    out = {}
    for name, arr in [
        ("rho", rho),
        ("u", u),
        ("v", v),
        ("p", p),
        ("eint", eint),
    ]:
        xu, au = compress_duplicate_x(x, arr)
        out[f"x_{name}"] = xu
        out[name] = au
    return out


def parse_exact(path: Path):
    x = []
    rho = []
    u = []
    p = []
    eint = []
    with path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            x.append(float(row["x"]))
            rho.append(float(row["rho"]))
            u.append(float(row["u"]))
            p.append(float(row["p"]))
            eint.append(float(row["eint"]))
    out = {}
    for name, arr in [("rho", rho), ("u", u), ("p", p), ("eint", eint)]:
        xa, aa = compress_duplicate_x(np.array(x), np.array(arr))
        out[f"x_{name}"] = xa
        out[name] = aa
    return out


def compare_profiles(xa, ya, xb, yb, n=4001):
    xlo = max(float(np.min(xa)), float(np.min(xb)))
    xhi = min(float(np.max(xa)), float(np.max(xb)))
    if not (xhi > xlo):
        raise RuntimeError("No overlap for profile comparison")
    xq = np.linspace(xlo, xhi, n)
    yaq = np.interp(xq, xa, ya)
    ybq = np.interp(xq, xb, yb)
    d = yaq - ybq
    return {
        "L1": float(np.mean(np.abs(d))),
        "L2": float(math.sqrt(np.mean(d * d))),
        "Linf": float(np.max(np.abs(d))),
        "xlo": xlo,
        "xhi": xhi,
    }


def run_fextract(fextract: Path, plotfile: Path, out_csv: Path, direction: int, x=None, y=None):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(fextract),
        "-d",
        str(direction),
        "-v",
        "rho momx momy E",
        "-csv",
        "-s",
        str(out_csv),
    ]
    if x is not None:
        cmd += ["-x", f"{x:.16g}"]
    if y is not None:
        cmd += ["-y", f"{y:.16g}"]
    cmd += [str(plotfile)]
    subprocess.run(cmd, check=True)


def draw_case_figure(case, out_dir: Path, tmp_dir: Path, rho_range, fsnapshot: Path, fboxinfo: Path):
    plotfile = case["plotfile"]
    geom = parse_plot_geom(plotfile)
    xlo, ylo = geom["prob_lo"][0], geom["prob_lo"][1]
    xhi, yhi = geom["prob_hi"][0], geom["prob_hi"][1]
    nx0, ny0 = geom["ncell"][0], geom["ncell"][1]
    test = case["test"]
    ic = case["ic"]
    mode = case["mode"]
    split = IC_LABEL[ic]
    rmin, rmax = rho_range

    out_ppm = tmp_dir / f"{plotfile.name}.rho.ppm"
    make_snapshot(plotfile, rmin, rmax, out_ppm, fsnapshot)
    img = read_ppm_p6(out_ppm)

    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "font.size": 10.5,
            "axes.titlesize": 11.5,
            "axes.labelsize": 10.5,
        }
    )

    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    ax.imshow(img, extent=(xlo, xhi, ylo, yhi), origin="upper", interpolation="nearest")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")

    mode_title = "Uniform (400x400)" if mode == "uniform" else "AMR (base 200x200, max_level=1)"
    ax.set_title(
        f"Toro test {test}, {IC_TITLE[ic]}, {mode_title}\n"
        f"rho, t = {geom['time']:.3g}, fixed scale [{rmin:.3f}, {rmax:.3f}]"
    )

    if mode == "amr":
        boxes = parse_boxes(plotfile, fboxinfo)
        dx0 = (xhi - xlo) / nx0
        dy0 = (yhi - ylo) / ny0
        for lev in sorted(boxes.keys()):
            if lev < 1:
                continue
            rr = 2**lev
            dx = dx0 / rr
            dy = dy0 / rr
            for ilo, jlo, ihi, jhi in boxes[lev]:
                xa = xlo + ilo * dx
                xb = xlo + (ihi + 1) * dx
                ya = ylo + jlo * dy
                yb = ylo + (jhi + 1) * dy
                ax.add_patch(
                    Rectangle((xa, ya), xb - xa, yb - ya, fill=False, edgecolor="white", linewidth=0.35, alpha=0.95)
                )

    out_dir.mkdir(parents=True, exist_ok=True)
    if mode == "uniform":
        stem = out_dir / f"rho_t{test}_{split}_uniform"
    else:
        stem = out_dir / f"rho_t{test}_{split}_amr_with_patches"
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=500, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return [stem.with_suffix(".png"), stem.with_suffix(".pdf")]


def write_metrics_md(rows, out_path: Path):
    rows_by_test = {}
    for r in rows:
        rows_by_test.setdefault(r["test"], []).append(r)

    lines = []
    lines.append("# Task 2 Compact Metrics by Test")
    lines.append("")
    lines.append("Metrics are computed from centerline profiles with interpolation to a common coordinate grid.")
    lines.append("Rotation check compares x-split against y-split (normal velocity used for y-split).")
    lines.append("")
    for test in sorted(rows_by_test.keys()):
        lines.append(f"## Toro test {test}")
        lines.append("")
        lines.append("| mode | metric | case | variable | L1 | L2 | Linf |")
        lines.append("|---|---|---|---:|---:|---:|---:|")
        for r in sorted(rows_by_test[test], key=lambda x: (x["mode"], x["metric"], x["case"], x["variable"])):
            lines.append(
                f"| {r['mode']} | {r['metric']} | {r['case']} | {r['variable']} | "
                f"{r['L1']:.3e} | {r['L2']:.3e} | {r['Linf']:.3e} |"
            )
        lines.append("")
    out_path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Analyze and plot full Task-2 2D matrix.")
    ap.add_argument("--results-dir", required=True, help="Directory produced by run_task2_2d_fullmatrix_report.sh")
    args = ap.parse_args()

    results_dir = Path(args.results_dir).resolve()
    checks_dir = results_dir / "checks"
    figs_main = results_dir / "figures" / "main_text"
    figs_appendix = results_dir / "figures" / "appendix"
    tmp_ppm_dir = checks_dir / "_tmp_ppm"
    slice_dir = checks_dir / "slices"
    for d in [checks_dir, figs_main, figs_appendix, tmp_ppm_dir, slice_dir]:
        d.mkdir(parents=True, exist_ok=True)

    fsnapshot = find_tool("fsnapshot")
    fextrema = find_tool("fextrema")
    fboxinfo = find_tool("fboxinfo")
    fextract = find_tool("fextract")

    run_index = parse_run_index(results_dir / "logs" / "run_index.tsv")
    run_map = {(r["test"], r["ic"], r["mode"]): r for r in run_index}

    # Per-test rho ranges for fair color normalization within each test.
    rho_ranges = {}
    with (checks_dir / "rho_extrema_by_test.tsv").open("w") as f:
        f.write("test\tmode\tic\tplotfile\trho_min\trho_max\n")
        for test in range(1, 6):
            rmin = float("inf")
            rmax = float("-inf")
            for ic in [0, 1, 2]:
                for mode in ["uniform", "amr"]:
                    case = run_map[(test, ic, mode)]
                    mn, mx = parse_fextrema_rho(case["plotfile"], fextrema)
                    f.write(f"{test}\t{mode}\t{ic}\t{case['plotfile'].name}\t{mn:.16g}\t{mx:.16g}\n")
                    rmin = min(rmin, mn)
                    rmax = max(rmax, mx)
            rho_ranges[test] = (rmin, rmax)
            f.write(f"{test}\tall\tall\tall\t{rmin:.16g}\t{rmax:.16g}\n")

    # Generate figures with body/appendix split.
    outputs = []
    for test in range(1, 6):
        out_dir = figs_main / "test1" if test == 1 else figs_appendix / f"test{test}"
        for ic in [0, 1, 2]:
            for mode in ["uniform", "amr"]:
                outputs.extend(
                    draw_case_figure(
                        run_map[(test, ic, mode)],
                        out_dir,
                        tmp_ppm_dir,
                        rho_ranges[test],
                        fsnapshot,
                        fboxinfo,
                    )
                )

    with (results_dir / "figures" / "figure_manifest_fullmatrix.txt").open("w") as f:
        for p in outputs:
            f.write(str(p) + "\n")

    # Extract centerline data for metric tables.
    for test in range(1, 6):
        for mode in ["uniform", "amr"]:
            run_fextract(
                fextract,
                run_map[(test, 0, mode)]["plotfile"],
                slice_dir / f"t{test}_ic0_{mode}_centerline.csv",
                direction=0,
                y=0.5,
            )
            run_fextract(
                fextract,
                run_map[(test, 1, mode)]["plotfile"],
                slice_dir / f"t{test}_ic1_{mode}_centerline.csv",
                direction=1,
                x=0.5,
            )
            run_fextract(
                fextract,
                run_map[(test, 2, mode)]["plotfile"],
                slice_dir / f"t{test}_ic2_{mode}_xslice.csv",
                direction=0,
                y=0.5,
            )
            run_fextract(
                fextract,
                run_map[(test, 2, mode)]["plotfile"],
                slice_dir / f"t{test}_ic2_{mode}_yslice.csv",
                direction=1,
                x=0.5,
            )

    exec_dir = Path(__file__).resolve().parents[1]
    metrics_rows = []
    for test in range(1, 6):
        exact = parse_exact(exec_dir / f"exact_toro2nd_t{test}_n400.csv")
        for mode in ["uniform", "amr"]:
            xprof = parse_slice(slice_dir / f"t{test}_ic0_{mode}_centerline.csv")
            yprof = parse_slice(slice_dir / f"t{test}_ic1_{mode}_centerline.csv")
            dpx = parse_slice(slice_dir / f"t{test}_ic2_{mode}_xslice.csv")
            dpy = parse_slice(slice_dir / f"t{test}_ic2_{mode}_yslice.csv")

            # Rotation-equivalence: compare x-split against y-split with normal velocity mapping.
            for var, yvar in [("rho", "rho"), ("u", "v"), ("p", "p"), ("eint", "eint")]:
                c = compare_profiles(xprof[f"x_{var}"], xprof[var], yprof[f"x_{yvar}"], yprof[yvar])
                metrics_rows.append(
                    {
                        "test": test,
                        "mode": mode,
                        "metric": "rotation_xy",
                        "case": "xsplit_vs_ysplit",
                        "variable": var,
                        **c,
                    }
                )

            # 1D exact consistency on centerline.
            for split_case, prof, vel_var in [("xsplit", xprof, "u"), ("ysplit", yprof, "v")]:
                for var, pvar in [("rho", "rho"), ("u", vel_var), ("p", "p"), ("eint", "eint")]:
                    c = compare_profiles(prof[f"x_{pvar}"], prof[pvar], exact[f"x_{var}"], exact[var])
                    metrics_rows.append(
                        {
                            "test": test,
                            "mode": mode,
                            "metric": "exact_1d",
                            "case": f"{split_case}_vs_exact",
                            "variable": var,
                            **c,
                        }
                    )

            # Diagonal sanity check in rho.
            c = compare_profiles(dpx["x_rho"], dpx["rho"], dpy["x_rho"], dpy["rho"])
            metrics_rows.append(
                {
                    "test": test,
                    "mode": mode,
                    "metric": "diag_xy",
                    "case": "diag_xslice_vs_yslice",
                    "variable": "rho",
                    **c,
                }
            )

    with (checks_dir / "task2_metrics_long.csv").open("w", newline="") as f:
        fieldnames = ["test", "mode", "metric", "case", "variable", "L1", "L2", "Linf", "xlo", "xhi"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in metrics_rows:
            w.writerow(r)

    write_metrics_md(metrics_rows, checks_dir / "task2_metrics_by_test.md")

    report_lines = []
    report_lines.append("# Task 2 Reporting Notes")
    report_lines.append("")
    report_lines.append("## Resolution Matching Statement")
    report_lines.append("- Uniform runs: `amr.n_cell = 400 400`, `amr.max_level = 0`.")
    report_lines.append("- AMR runs: `amr.n_cell = 200 200`, `amr.max_level = 1`, `amr.ref_ratio = 2`.")
    report_lines.append("- Effective finest-grid spacing is matched (`~1/400`) between uniform and AMR.")
    report_lines.append("- `amr.blocking_factor = 8` divides both 400 and 200 exactly.")
    report_lines.append("")
    report_lines.append("## Figure Organization for Write-Up")
    report_lines.append("- Main text figures: `figures/main_text/test1/`.")
    report_lines.append("- Appendix figures (additional tests): `figures/appendix/test2/` to `figures/appendix/test5/`.")
    report_lines.append("- Within each test, six figures are provided: `ic=0/1/2`, each with uniform and AMR.")
    report_lines.append("")
    report_lines.append("## Metric Files")
    report_lines.append("- Compact markdown tables: `checks/task2_metrics_by_test.md`.")
    report_lines.append("- Full CSV metrics: `checks/task2_metrics_long.csv`.")
    report_lines.append("- Per-test rho normalization: `checks/rho_extrema_by_test.tsv`.")
    (checks_dir / "task2_report_notes.md").write_text("\n".join(report_lines) + "\n")

    print(f"RESULTS_DIR={results_dir}")
    print(f"FIGURE_MANIFEST={results_dir / 'figures' / 'figure_manifest_fullmatrix.txt'}")
    print(f"METRICS_MD={checks_dir / 'task2_metrics_by_test.md'}")
    print(f"METRICS_CSV={checks_dir / 'task2_metrics_long.csv'}")


if __name__ == "__main__":
    main()
