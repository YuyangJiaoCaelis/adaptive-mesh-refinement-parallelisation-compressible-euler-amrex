// Project-specific problem setup for the public project version
#include "Prob.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

namespace {
// Component indices
constexpr int RHO = 0;
constexpr int MX  = 1;
constexpr int MY  = 2;
constexpr amrex::Real TWO_PI = 6.2831853071795864769;
#if (AMREX_SPACEDIM == 3)
constexpr int MZ  = 3;
constexpr int E   = 4;
#else
constexpr int E   = 3;
#endif

struct TestState {
    amrex::Real rhoL, uL, vL, pL;
    amrex::Real rhoR, uR, vR, pR;
};

enum class InitMode : int {
    XSplit   = 0,  // 1D-style jump at x = x0
    YSplit   = 1,  // 1D-style jump at y = y0
    Diagonal = 2,  // 1D-style jump across n.(x-x0,y-y0) = 0
    Quadrant = 3,  // Lax-Liu 2D quadrant problem
    SmoothEntropy = 4 // periodic smooth entropy-wave advection
};

enum class SmoothProfile : int {
    Sinusoid = 0,
    Gaussian = 1
};

TestState get_test_state(int test_id)
{
    // Toro Table 4.1, gamma = 1.4
    switch (test_id) {
        case 1:
            return {1.0, 0.0, 0.0, 1.0, 0.125, 0.0, 0.0, 0.1};
        case 2:
            return {1.0, -2.0, 0.0, 0.4, 1.0,  2.0, 0.0, 0.4};
        case 3:
            return {1.0, 0.0, 0.0, 1000.0, 1.0, 0.0, 0.0, 0.01};
        case 4:
            return {1.0, 0.0, 0.0, 0.01, 1.0, 0.0, 0.0, 100.0};
        case 5:
            return {5.99924, 19.5975, 0.0, 460.894, 5.99242, -6.19633, 0.0, 46.0950};
        default:
            return {1.0, 0.0, 0.0, 1.0, 0.125, 0.0, 0.0, 0.1};
    }
}

} // namespace

void initdata(amrex::MultiFab& state, amrex::Geometry const& geom, amrex::Real /*time*/)
{
    static bool inited = false;
    static int test_id = 1;
    static int ic_mode = static_cast<int>(InitMode::XSplit);
    static int print_ic = 0;
    static amrex::Real gamma = 1.4;
    static amrex::Real x0 = 0.5;
    static amrex::Real y0 = 0.5;
    static amrex::Real diag_nx = 1.0;
    static amrex::Real diag_ny = 1.0;
    static amrex::Real smooth_rho0 = 1.0;
    static amrex::Real smooth_amp = 0.2;
    static amrex::Real smooth_u = 1.0;
    static amrex::Real smooth_v = 1.0;
    static amrex::Real smooth_p = 1.0;
    static amrex::Real smooth_kx = 1.0;
    static amrex::Real smooth_ky = 1.0;
    static int smooth_profile = static_cast<int>(SmoothProfile::Sinusoid);
    static amrex::Real smooth_sigma = 0.08;

    if (!inited) {
        amrex::ParmParse pp("prob");
        pp.query("test", test_id);
        pp.query("ic", ic_mode);
        pp.query("print_ic", print_ic);
        pp.query("gamma", gamma);
        pp.query("x0", x0);
        pp.query("y0", y0);
        pp.query("diag_nx", diag_nx);
        pp.query("diag_ny", diag_ny);
        pp.query("smooth_rho0", smooth_rho0);
        pp.query("smooth_amp", smooth_amp);
        pp.query("smooth_u", smooth_u);
        pp.query("smooth_v", smooth_v);
        pp.query("smooth_p", smooth_p);
        pp.query("smooth_kx", smooth_kx);
        pp.query("smooth_ky", smooth_ky);
        pp.query("smooth_profile", smooth_profile);
        pp.query("smooth_sigma", smooth_sigma);

        const amrex::Real nrm = std::sqrt(diag_nx*diag_nx + diag_ny*diag_ny);
        if (nrm > 0.0) {
            diag_nx /= nrm;
            diag_ny /= nrm;
        } else {
            diag_nx = 1.0;
            diag_ny = 0.0;
        }

        if (print_ic != 0) {
            const TestState ts_print = get_test_state(test_id);
            amrex::Print()
                << "Initial condition: test=" << test_id
                << " prob.ic=" << ic_mode
                << " gamma=" << gamma
                << " | L: (rho,u,v,p)=("
                << ts_print.rhoL << ", " << ts_print.uL << ", " << ts_print.vL << ", " << ts_print.pL << ")"
                << " R: (rho,u,v,p)=("
                << ts_print.rhoR << ", " << ts_print.uR << ", " << ts_print.vR << ", " << ts_print.pR << ")"
                << "\n";
            if (ic_mode == static_cast<int>(InitMode::Quadrant)) {
                amrex::Print() << "Initial condition: Lax-Liu quadrant states\n";
            } else if (ic_mode == static_cast<int>(InitMode::SmoothEntropy)) {
                amrex::Print()
                    << "Initial condition: smooth entropy-wave data"
                    << " rho0=" << smooth_rho0
                    << " amp=" << smooth_amp
                    << " u=" << smooth_u
                    << " v=" << smooth_v
                    << " p=" << smooth_p
                    << " profile=" << smooth_profile
                    << " k=(" << smooth_kx << "," << smooth_ky << ")"
                    << " sigma=" << smooth_sigma << "\n";
            }
        }

        inited = true;
    }

    const auto prob_lo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

    const TestState ts = get_test_state(test_id);

    for (amrex::MFIter mfi(state); mfi.isValid(); ++mfi) {
        auto const& box = mfi.validbox();
        auto const& arr = state.array(mfi);

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            amrex::Real rho, u, v, p;

            if (ic_mode == static_cast<int>(InitMode::SmoothEntropy)) {
                if (smooth_profile == static_cast<int>(SmoothProfile::Gaussian)) {
                    const amrex::Real dxg = x - x0;
                    const amrex::Real dyg = y - y0;
                    const amrex::Real r2 = dxg*dxg + dyg*dyg;
                    const amrex::Real inv_sigma2 = 1.0 / (smooth_sigma * smooth_sigma);
                    rho = smooth_rho0 + smooth_amp * std::exp(-r2 * inv_sigma2);
                } else {
                    const amrex::Real phase =
                        TWO_PI * (smooth_kx * (x - x0) + smooth_ky * (y - y0));
                    rho = smooth_rho0 + smooth_amp * std::sin(phase);
                }
                u = smooth_u;
                v = smooth_v;
                p = smooth_p;
            } else if (ic_mode == static_cast<int>(InitMode::Quadrant)) {
                // Lax & Liu (1998) quadrant problem, split by x=x0 and y=y0.
                if (x >= x0 && y >= y0) {
                    rho = 1.5;    u = 0.0;    v = 0.0;    p = 1.5;
                } else if (x < x0 && y >= y0) {
                    rho = 0.5323; u = 1.206;  v = 0.0;    p = 0.3;
                } else if (x < x0 && y < y0) {
                    rho = 0.138;  u = 1.206;  v = 1.206;  p = 0.029;
                } else {
                    rho = 0.5323; u = 0.0;    v = 1.206;  p = 0.3;
                }
            } else {
                // Treat Toro states as (rho, un, ut, p) where un is velocity
                // normal to the initial discontinuity and ut is tangential.
                amrex::Real nx = 1.0;
                amrex::Real ny = 0.0;
                amrex::Real signed_dist = x - x0;

                if (ic_mode == static_cast<int>(InitMode::YSplit)) {
                    nx = 0.0;
                    ny = 1.0;
                    signed_dist = y - y0;
                } else if (ic_mode == static_cast<int>(InitMode::Diagonal)) {
                    nx = diag_nx;
                    ny = diag_ny;
                    signed_dist = nx*(x - x0) + ny*(y - y0);
                }

                const bool left_state = (signed_dist < 0.0);
                amrex::Real un, ut;
                if (left_state) {
                    rho = ts.rhoL; un = ts.uL; ut = ts.vL; p = ts.pL;
                } else {
                    rho = ts.rhoR; un = ts.uR; ut = ts.vR; p = ts.pR;
                }

                // n=(nx,ny), t=(-ny,nx): (u,v) = un*n + ut*t.
                u = un*nx - ut*ny;
                v = un*ny + ut*nx;
            }

            arr(i,j,k,RHO) = rho;
            arr(i,j,k,MX)  = rho * u;
            arr(i,j,k,MY)  = rho * v;
#if (AMREX_SPACEDIM == 3)
            arr(i,j,k,MZ)  = 0.0;
#endif
            const amrex::Real kinetic = 0.5 * rho * (u*u + v*v);
            arr(i,j,k,E)   = p/(gamma - 1.0) + kinetic;
        });
    }
}
