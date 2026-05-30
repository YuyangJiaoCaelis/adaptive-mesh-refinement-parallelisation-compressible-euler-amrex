#!/usr/bin/env python3
import argparse
import math
import os
import re
import subprocess
from pathlib import Path


def parse_inputs(path):
    prob = {"test": 1, "gamma": 1.4, "x0": 0.5, "stop_time": None}
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # strip inline comments
            if "#" in line:
                line = line.split("#", 1)[0].strip()
            if line.startswith("prob.test"):
                prob["test"] = int(line.split("=", 1)[1].strip())
            elif line.startswith("prob.gamma"):
                prob["gamma"] = float(line.split("=", 1)[1].strip())
            elif line.startswith("prob.x0"):
                prob["x0"] = float(line.split("=", 1)[1].strip())
            elif line.startswith("stop_time"):
                prob["stop_time"] = float(line.split("=", 1)[1].strip())
    return prob


def run_fextract(plotfile, slice_path):
    if "AMREX_HOME" not in os.environ:
        raise EnvironmentError("AMREX_HOME must be set so fextract.gnu.ex can be found")
    amrex_home = Path(os.environ["AMREX_HOME"])
    fextract = amrex_home / "Tools" / "Plotfile" / "fextract.gnu.ex"
    if not fextract.exists():
        raise FileNotFoundError(f"fextract not found at {fextract}")
    cmd = [
        str(fextract),
        "-s", str(slice_path),
        "-v", "rho momx momy E",
        str(plotfile),
    ]
    subprocess.run(cmd, check=True)


def read_slice(slice_path):
    t = None
    rows = []
    with open(slice_path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("# time"):
                    # format: "# time = 0.25"
                    parts = line.split("=")
                    if len(parts) >= 2:
                        try:
                            t = float(parts[1].strip())
                        except ValueError:
                            pass
                continue
            parts = line.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) < 5:
                continue
            rows.append(vals[:5])  # x, rho, momx, momy, E
    return t, rows


def infer_plotfile_from_slice(slice_path):
    sp = Path(slice_path)
    if sp.suffix == ".slice":
        candidate = Path(str(sp)[:-len(".slice")])
        if candidate.is_dir():
            return candidate
    return None


def parse_plotfile_header(plotfile):
    header_path = Path(plotfile) / "Header"
    if not header_path.exists():
        raise FileNotFoundError(f"Plotfile header not found at {header_path}")

    with open(header_path) as f:
        lines = [ln.strip() for ln in f if ln.strip() != ""]

    if len(lines) < 8:
        raise SystemExit(f"Header too short: {header_path}")

    idx = 0
    _version = lines[idx]
    idx += 1
    nvars = int(lines[idx])
    idx += 1
    var_names = lines[idx:idx + nvars]
    idx += nvars
    dim = int(lines[idx])
    idx += 1
    time = float(lines[idx])
    idx += 1
    finest_level = int(lines[idx])
    idx += 1
    prob_lo = [float(x) for x in lines[idx].split()]
    idx += 1
    prob_hi = [float(x) for x in lines[idx].split()]
    idx += 1

    if len(prob_lo) != dim or len(prob_hi) != dim:
        raise SystemExit(f"Unexpected prob_lo/prob_hi size in {header_path}")

    # Read level-0 domain line: e.g. "((0,0) (199,3) (0,0))"
    if idx >= len(lines):
        raise SystemExit(f"Missing level domain in {header_path}")
    dom_line = lines[idx]
    idx += 1
    nums = [int(v) for v in re.findall(r"-?\d+", dom_line)]
    if len(nums) < 2 * dim:
        raise SystemExit(f"Could not parse level domain from line: {dom_line}")
    dom_lo = nums[0:dim]
    dom_hi = nums[dim:2 * dim]
    ncell = [dom_hi[d] - dom_lo[d] + 1 for d in range(dim)]

    # Find first dim-length positive-float line after domain: this is level-0 dx.
    dx = None
    for line in lines[idx:]:
        parts = line.split()
        if len(parts) != dim:
            continue
        ok = True
        vals = []
        for p in parts:
            try:
                vals.append(float(p))
            except ValueError:
                ok = False
                break
        if ok and all(v > 0.0 for v in vals):
            dx = vals
            break
    if dx is None:
        raise SystemExit(f"Could not find level-0 dx in {header_path}")

    return {
        "header_path": str(header_path),
        "var_names": var_names,
        "dim": dim,
        "time": time,
        "finest_level": finest_level,
        "prob_lo": prob_lo,
        "prob_hi": prob_hi,
        "dom_lo": dom_lo,
        "dom_hi": dom_hi,
        "ncell": ncell,
        "dx": dx,
    }


def check_x_coords_against_geometry(xs, geom, atol=1.0e-12):
    if geom["dim"] < 1:
        raise SystemExit("Invalid geometry dimension")
    nx = geom["ncell"][0]
    if len(xs) != nx:
        raise SystemExit(f"x-size mismatch: slice has {len(xs)} points, header nx={nx}")

    x0 = geom["prob_lo"][0]
    dx = geom["dx"][0]
    max_abs = 0.0
    for i, xv in enumerate(xs):
        xc = x0 + (i + 0.5) * dx
        max_abs = max(max_abs, abs(xv - xc))
    if max_abs > atol:
        raise SystemExit(f"x-centers mismatch against header geometry: max |dx|={max_abs}")
    return max_abs


def write_derived(out_path, rows, gamma):
    with open(out_path, "w") as f:
        f.write("x,rho,u,v,p,eint\n")
        for x, rho, momx, momy, E in rows:
            if rho <= 0.0:
                u = 0.0
                v = 0.0
            else:
                u = momx / rho
                v = momy / rho
            kinetic = 0.5 * rho * (u * u + v * v)
            p = (gamma - 1.0) * (E - kinetic)
            if rho > 0.0:
                eint = p / ((gamma - 1.0) * rho)
            else:
                eint = 0.0
            f.write(f"{x},{rho},{u},{v},{p},{eint}\n")


def read_derived(path):
    rows = []
    with open(path) as f:
        header = f.readline().strip().split(",")
        if header[:6] != ["x", "rho", "u", "v", "p", "eint"]:
            raise SystemExit(f"Unexpected header in {path}: {header}")
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 6:
                continue
            rows.append([float(x) for x in parts[:6]])
    return rows


def verify_derived(out_path, rows, gamma, rtol=1.0e-9, atol=1.0e-9):
    derived = read_derived(out_path)
    if len(derived) != len(rows):
        raise SystemExit(f"Derived rows mismatch: {len(derived)} vs {len(rows)}")

    max_du = 0.0
    max_dp = 0.0
    for (x, rho, momx, momy, E), (x2, rho2, u, v, p, _eint) in zip(rows, derived):
        if abs(x - x2) > 1.0e-12:
            raise SystemExit(f"x mismatch at {x} vs {x2}")
        if abs(rho - rho2) > 1.0e-12:
            raise SystemExit(f"rho mismatch at x={x}: {rho} vs {rho2}")
        if rho <= 0.0:
            u_chk = 0.0
        else:
            u_chk = momx / rho
        p_chk = (gamma - 1.0) * (E - 0.5 * (momx * momx + momy * momy) / rho)

        du = abs(u - u_chk)
        dp = abs(p - p_chk)
        max_du = max(max_du, du)
        max_dp = max(max_dp, dp)

        if du > max(atol, rtol * max(1.0, abs(u_chk))):
            raise SystemExit(f"u mismatch at x={x}: u={u} vs u_chk={u_chk}")
        if dp > max(atol, rtol * max(1.0, abs(p_chk))):
            raise SystemExit(f"p mismatch at x={x}: p={p} vs p_chk={p_chk}")

    print(f"Derived check passed: max |du|={max_du:.3e}, max |dp|={max_dp:.3e}")


# Toro exact Riemann solver (ideal gas)
TESTS = {
    1: dict(rhoL=1.0, uL=0.0, pL=1.0,   rhoR=0.125, uR=0.0,  pR=0.1),
    2: dict(rhoL=1.0, uL=-2.0, pL=0.4,  rhoR=1.0,   uR=2.0,  pR=0.4),
    3: dict(rhoL=1.0, uL=0.0, pL=1000.0, rhoR=1.0,  uR=0.0,  pR=0.01),
    4: dict(rhoL=1.0, uL=0.0, pL=0.01, rhoR=1.0,   uR=0.0,  pR=100.0),
    5: dict(rhoL=5.99924, uL=19.5975, pL=460.894, rhoR=5.99242, uR=-6.19633, pR=46.0950),
}


def sound_speed(gamma, p, rho):
    return math.sqrt(gamma * p / rho)


def f_pressure(p, rhoK, uK, pK, aK, gamma):
    if p <= pK:
        pr = p / pK
        f = (2 * aK / (gamma - 1.0)) * (pr ** ((gamma - 1.0) / (2 * gamma)) - 1.0)
        df = (1.0 / (rhoK * aK)) * pr ** (-(gamma + 1.0) / (2 * gamma))
    else:
        A = 2.0 / ((gamma + 1.0) * rhoK)
        B = (gamma - 1.0) / (gamma + 1.0) * pK
        sqrt_term = math.sqrt(A / (p + B))
        f = (p - pK) * sqrt_term
        df = sqrt_term * (1.0 - 0.5 * (p - pK) / (p + B))
    return f, df


def solve_star_pressure(rhoL, uL, pL, rhoR, uR, pR, gamma):
    aL = sound_speed(gamma, pL, rhoL)
    aR = sound_speed(gamma, pR, rhoR)
    pPV = 0.5 * (pL + pR) - 0.125 * (uR - uL) * (rhoL + rhoR) * (aL + aR)
    p = max(1.0e-8, pPV)
    for _ in range(50):
        fL, dfL = f_pressure(p, rhoL, uL, pL, aL, gamma)
        fR, dfR = f_pressure(p, rhoR, uR, pR, aR, gamma)
        f = fL + fR + (uR - uL)
        df = dfL + dfR
        p_new = p - f / df
        if p_new < 0:
            p_new = 1.0e-8
        if abs(p_new - p) / (p_new + 1.0e-12) < 1.0e-8:
            p = p_new
            break
        p = p_new
    return p


def star_velocity(p_star, rhoL, uL, pL, rhoR, uR, pR, gamma):
    aL = sound_speed(gamma, pL, rhoL)
    aR = sound_speed(gamma, pR, rhoR)
    fL, _ = f_pressure(p_star, rhoL, uL, pL, aL, gamma)
    fR, _ = f_pressure(p_star, rhoR, uR, pR, aR, gamma)
    return 0.5 * (uL + uR) + 0.5 * (fR - fL)


def sample_solution(x, t, x0, left, right, gamma, p_star, u_star):
    rhoL, uL, pL = left
    rhoR, uR, pR = right
    aL = sound_speed(gamma, pL, rhoL)
    aR = sound_speed(gamma, pR, rhoR)

    if t == 0.0:
        return (rhoL, uL, pL) if x < x0 else (rhoR, uR, pR)

    xi = (x - x0) / t

    if xi <= u_star:
        if p_star <= pL:
            a_star_L = aL * (p_star / pL) ** ((gamma - 1.0) / (2 * gamma))
            SHL = uL - aL
            STL = u_star - a_star_L
            if xi < SHL:
                return rhoL, uL, pL
            elif xi > STL:
                rho_star_L = rhoL * (p_star / pL) ** (1.0 / gamma)
                return rho_star_L, u_star, p_star
            else:
                u = (2.0 / (gamma + 1.0)) * (aL + 0.5 * (gamma - 1.0) * uL + xi)
                a = (2.0 / (gamma + 1.0)) * (aL + 0.5 * (gamma - 1.0) * (uL - xi))
                rho = rhoL * (a / aL) ** (2.0 / (gamma - 1.0))
                p = pL * (a / aL) ** (2.0 * gamma / (gamma - 1.0))
                return rho, u, p
        else:
            SL = uL - aL * math.sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / pL - 1.0))
            if xi < SL:
                return rhoL, uL, pL
            rho_star_L = rhoL * ((p_star / pL + (gamma - 1.0) / (gamma + 1.0)) /
                                 ((gamma - 1.0) / (gamma + 1.0) * p_star / pL + 1.0))
            return rho_star_L, u_star, p_star
    else:
        if p_star <= pR:
            a_star_R = aR * (p_star / pR) ** ((gamma - 1.0) / (2 * gamma))
            SHR = uR + aR
            STR = u_star + a_star_R
            if xi > SHR:
                return rhoR, uR, pR
            elif xi < STR:
                rho_star_R = rhoR * (p_star / pR) ** (1.0 / gamma)
                return rho_star_R, u_star, p_star
            else:
                u = (2.0 / (gamma + 1.0)) * (-aR + 0.5 * (gamma - 1.0) * uR + xi)
                a = (2.0 / (gamma + 1.0)) * (aR - 0.5 * (gamma - 1.0) * (uR - xi))
                rho = rhoR * (a / aR) ** (2.0 / (gamma - 1.0))
                p = pR * (a / aR) ** (2.0 * gamma / (gamma - 1.0))
                return rho, u, p
        else:
            SR = uR + aR * math.sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / pR - 1.0))
            if xi > SR:
                return rhoR, uR, pR
            rho_star_R = rhoR * ((p_star / pR + (gamma - 1.0) / (gamma + 1.0)) /
                                 ((gamma - 1.0) / (gamma + 1.0) * p_star / pR + 1.0))
            return rho_star_R, u_star, p_star


def exact_metadata(test, gamma):
    if test not in TESTS:
        raise SystemExit(f"Unknown test {test}")
    st = TESTS[test]
    rhoL, uL, pL = st["rhoL"], st["uL"], st["pL"]
    rhoR, uR, pR = st["rhoR"], st["uR"], st["pR"]
    p_star = solve_star_pressure(rhoL, uL, pL, rhoR, uR, pR, gamma)
    u_star = star_velocity(p_star, rhoL, uL, pL, rhoR, uR, pR, gamma)
    return (rhoL, uL, pL), (rhoR, uR, pR), p_star, u_star


def wave_positions(t, x0, left, right, gamma, p_star, u_star):
    rhoL, uL, pL = left
    rhoR, uR, pR = right
    aL = sound_speed(gamma, pL, rhoL)
    aR = sound_speed(gamma, pR, rhoR)

    waves = []
    if p_star <= pL:
        a_star_L = aL * (p_star / pL) ** ((gamma - 1.0) / (2.0 * gamma))
        s_head = uL - aL
        s_tail = u_star - a_star_L
        waves.append(("left_rarefaction_head", s_head, x0 + s_head * t))
        waves.append(("left_rarefaction_tail", s_tail, x0 + s_tail * t))
    else:
        s_left = uL - aL * math.sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / pL - 1.0))
        waves.append(("left_shock", s_left, x0 + s_left * t))

    waves.append(("contact", u_star, x0 + u_star * t))

    if p_star <= pR:
        a_star_R = aR * (p_star / pR) ** ((gamma - 1.0) / (2.0 * gamma))
        s_tail = u_star + a_star_R
        s_head = uR + aR
        waves.append(("right_rarefaction_tail", s_tail, x0 + s_tail * t))
        waves.append(("right_rarefaction_head", s_head, x0 + s_head * t))
    else:
        s_right = uR + aR * math.sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (p_star / pR - 1.0))
        waves.append(("right_shock", s_right, x0 + s_right * t))

    return waves


def write_exact(out_path, xs, t, test, gamma, x0):
    left, right, p_star, u_star = exact_metadata(test, gamma)
    waves = wave_positions(t, x0, left, right, gamma, p_star, u_star)

    with open(out_path, "w") as f:
        f.write("x,rho,u,p,eint\n")
        for x in xs:
            rho, u, p = sample_solution(x, t, x0, left, right, gamma, p_star, u_star)
            eint = p / ((gamma - 1.0) * rho)
            f.write(f"{x},{rho},{u},{p},{eint}\n")
    return waves


def wave_alignment_from_rows(rows, waves):
    if len(rows) < 3:
        return []
    xs = [r[0] for r in rows]
    rho = [r[1] for r in rows]
    edges = []
    for i in range(len(xs) - 1):
        xedge = 0.5 * (xs[i] + xs[i + 1])
        drho = (rho[i + 1] - rho[i]) / (xs[i + 1] - xs[i])
        edges.append((xedge, abs(drho)))
    if not edges:
        return []

    out = []
    for name, speed, xpos in waves:
        xbest = min(edges, key=lambda e: abs(e[0] - xpos))[0]
        out.append((name, speed, xpos, xbest, abs(xbest - xpos)))
    return out


def main():
    parser = argparse.ArgumentParser(description="Postprocess 1D Euler plotfiles")
    parser.add_argument("--plotfile", type=str, help="Plotfile directory (pltXXXXX)")
    parser.add_argument("--slice", type=str, default=None, help="Optional precomputed slice file")
    parser.add_argument("--out", type=str, required=True, help="Output derived CSV")
    parser.add_argument("--inputs", type=str, default="inputs", help="Inputs file")
    parser.add_argument("--exact", type=str, default=None, help="Write exact solution CSV")
    parser.add_argument("--time", type=float, default=None, help="Override time for exact solution")
    args = parser.parse_args()

    base = Path(__file__).resolve().parent
    inputs_path = (base / args.inputs).resolve()
    prob = parse_inputs(inputs_path)

    plotfile_path = None
    if args.plotfile:
        plotfile_path = Path(args.plotfile)

    if args.slice:
        slice_path = Path(args.slice)
        if plotfile_path is None:
            plotfile_path = infer_plotfile_from_slice(slice_path)
    else:
        if plotfile_path is None:
            raise SystemExit("Provide --plotfile or --slice")
        slice_path = Path(str(plotfile_path) + ".slice")
        run_fextract(plotfile_path, slice_path)

    t_slice, rows = read_slice(slice_path)
    if not rows:
        raise SystemExit(f"No data rows found in {slice_path}")

    geom = None
    if plotfile_path is not None and plotfile_path.is_dir():
        geom = parse_plotfile_header(plotfile_path)
        x_err = check_x_coords_against_geometry([r[0] for r in rows], geom)
        print("Geometry check:"
              f" prob_lo={geom['prob_lo']} prob_hi={geom['prob_hi']}"
              f" ncell={geom['ncell']} dx={geom['dx']}"
              f" max|x_slice-x_cellcenter|={x_err:.3e}")

    out_path = Path(args.out)
    write_derived(out_path, rows, prob["gamma"])
    verify_derived(out_path, rows, prob["gamma"])

    if args.exact:
        t = args.time
        if t is None:
            if geom is not None:
                t = geom["time"]
            elif t_slice is not None:
                t = t_slice
            else:
                t = prob.get("stop_time", None)
        if t is None:
            raise SystemExit("Could not determine time for exact solution")
        xs = [r[0] for r in rows]
        waves = write_exact(args.exact, xs, t, prob["test"], prob["gamma"], prob["x0"])
        print(f"Exact setup: t={t} gamma={prob['gamma']} x0={prob['x0']} test={prob['test']}")
        print("Predicted wave locations (name, speed, x):")
        for name, speed, xpos in waves:
            print(f"  {name:24s} speed={speed: .12e} x={xpos: .12e}")
        aligns = wave_alignment_from_rows(rows, waves)
        print("Nearest numerical rho-gradient edges to predicted waves:")
        for name, speed, xpos, xnum, dxabs in aligns:
            _ = speed
            print(f"  {name:24s} x_exact={xpos: .12e} x_num~={xnum: .12e} |dx|={dxabs: .12e}")

    print(f"Wrote {out_path}")
    if args.exact:
        print(f"Wrote {args.exact}")


if __name__ == "__main__":
    main()
