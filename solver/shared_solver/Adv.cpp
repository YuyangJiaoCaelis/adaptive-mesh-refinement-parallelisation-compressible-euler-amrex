// Modified for the Euler AMR project - public project version
#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Math.H>
#include <AMReX_MFParallelFor.H>
#include <AMReX_Reduce.H>

#include "AmrLevelAdv.H"

#if (AMREX_SPACEDIM != 2)
#error "Euler_AmrLevel currently implemented for 2D only"
#endif

namespace {
constexpr int RHO = 0;
constexpr int MX  = 1;
constexpr int MY  = 2;
constexpr int E   = 3;
constexpr int NQ  = 4; // rho, u, v, p
constexpr amrex::Real RHO_FLOOR = 1.0e-12;
constexpr amrex::Real P_FLOOR   = 1.0e-12;

struct ToroState {
    amrex::Real rhoL;
    amrex::Real uL;
    amrex::Real vL;
    amrex::Real pL;
    amrex::Real rhoR;
    amrex::Real uR;
    amrex::Real vR;
    amrex::Real pR;
};

AMREX_GPU_HOST AMREX_FORCE_INLINE
ToroState get_toro_state (int test_id)
{
    switch (test_id) {
        case 1: return {1.0, 0.0, 0.0, 1.0, 0.125, 0.0, 0.0, 0.1};
        case 2: return {1.0, -2.0, 0.0, 0.4, 1.0, 2.0, 0.0, 0.4};
        case 3: return {1.0, 0.0, 0.0, 1000.0, 1.0, 0.0, 0.0, 0.01};
        case 4: return {1.0, 0.0, 0.0, 0.01, 1.0, 0.0, 0.0, 100.0};
        case 5: return {5.99924, 19.5975, 0.0, 460.894, 5.99242, -6.19633, 0.0, 46.0950};
        default: return {1.0, 0.0, 0.0, 1.0, 0.125, 0.0, 0.0, 0.1};
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real minmod(amrex::Real a, amrex::Real b)
{
    if (a*b <= 0.0) return 0.0;
    return (amrex::Math::abs(a) < amrex::Math::abs(b)) ? a : b;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void cons_to_prim(const amrex::Real* U, amrex::Real gamma,
                  amrex::Real& rho, amrex::Real& u, amrex::Real& v, amrex::Real& p)
{
    rho = amrex::max(U[RHO], RHO_FLOOR);
    u   = U[MX]/rho;
    v   = U[MY]/rho;
    const amrex::Real kinetic = 0.5*rho*(u*u + v*v);
    p   = (gamma-1.0)*(U[E] - kinetic);
    if (p < P_FLOOR) p = P_FLOOR;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void prim_to_cons(amrex::Real rho, amrex::Real u, amrex::Real v, amrex::Real p,
                  amrex::Real gamma, amrex::Real* U)
{
    U[RHO] = rho;
    U[MX]  = rho*u;
    U[MY]  = rho*v;
    U[E]   = p/(gamma-1.0) + 0.5*rho*(u*u + v*v);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void flux_x_from_prim(amrex::Real rho, amrex::Real u, amrex::Real v, amrex::Real p,
                      amrex::Real gamma, amrex::Real* F)
{
    const amrex::Real Etot = p/(gamma-1.0) + 0.5*rho*(u*u + v*v);
    F[RHO] = rho*u;
    F[MX]  = rho*u*u + p;
    F[MY]  = rho*u*v;
    F[E]   = u*(Etot + p);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void flux_y_from_prim(amrex::Real rho, amrex::Real u, amrex::Real v, amrex::Real p,
                      amrex::Real gamma, amrex::Real* G)
{
    const amrex::Real Etot = p/(gamma-1.0) + 0.5*rho*(u*u + v*v);
    G[RHO] = rho*v;
    G[MX]  = rho*u*v;
    G[MY]  = rho*v*v + p;
    G[E]   = v*(Etot + p);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hll_flux_x(amrex::Real rhoL, amrex::Real uL, amrex::Real vL, amrex::Real pL,
                amrex::Real rhoR, amrex::Real uR, amrex::Real vR, amrex::Real pR,
                amrex::Real gamma, amrex::Real* F)
{
    const amrex::Real cL = std::sqrt(gamma*pL/rhoL);
    const amrex::Real cR = std::sqrt(gamma*pR/rhoR);
    const amrex::Real SL = amrex::min(uL - cL, uR - cR);
    const amrex::Real SR = amrex::max(uL + cL, uR + cR);

    amrex::Real UL[4], UR[4], FL[4], FR[4];
    prim_to_cons(rhoL,uL,vL,pL,gamma,UL);
    prim_to_cons(rhoR,uR,vR,pR,gamma,UR);
    flux_x_from_prim(rhoL,uL,vL,pL,gamma,FL);
    flux_x_from_prim(rhoR,uR,vR,pR,gamma,FR);

    if (SL >= 0.0) {
        for (int n=0; n<4; ++n) F[n] = FL[n];
    } else if (SR <= 0.0) {
        for (int n=0; n<4; ++n) F[n] = FR[n];
    } else {
        const amrex::Real inv = 1.0/(SR - SL);
        for (int n=0; n<4; ++n) {
            F[n] = (SR*FL[n] - SL*FR[n] + SL*SR*(UR[n]-UL[n]))*inv;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void hll_flux_y(amrex::Real rhoL, amrex::Real uL, amrex::Real vL, amrex::Real pL,
                amrex::Real rhoR, amrex::Real uR, amrex::Real vR, amrex::Real pR,
                amrex::Real gamma, amrex::Real* G)
{
    const amrex::Real cL = std::sqrt(gamma*pL/rhoL);
    const amrex::Real cR = std::sqrt(gamma*pR/rhoR);
    const amrex::Real SL = amrex::min(vL - cL, vR - cR);
    const amrex::Real SR = amrex::max(vL + cL, vR + cR);

    amrex::Real UL[4], UR[4], GL[4], GR[4];
    prim_to_cons(rhoL,uL,vL,pL,gamma,UL);
    prim_to_cons(rhoR,uR,vR,pR,gamma,UR);
    flux_y_from_prim(rhoL,uL,vL,pL,gamma,GL);
    flux_y_from_prim(rhoR,uR,vR,pR,gamma,GR);

    if (SL >= 0.0) {
        for (int n=0; n<4; ++n) G[n] = GL[n];
    } else if (SR <= 0.0) {
        for (int n=0; n<4; ++n) G[n] = GR[n];
    } else {
        const amrex::Real inv = 1.0/(SR - SL);
        for (int n=0; n<4; ++n) {
            G[n] = (SR*GL[n] - SL*GR[n] + SL*SR*(UR[n]-UL[n]))*inv;
        }
    }
}

} // namespace

AMREX_GPU_HOST
void
AmrLevelAdv::advect (const amrex::Real time,
                     const amrex::Box& bx,
                     amrex::GpuArray<amrex::Box,BL_SPACEDIM> nbx,
                     const amrex::FArrayBox& statein,
                     amrex::FArrayBox& stateout,
                     AMREX_D_DECL(const amrex::FArrayBox& /*vx*/,
                                  const amrex::FArrayBox& /*vy*/,
                                  const amrex::FArrayBox& /*vz*/),
                     AMREX_D_DECL(amrex::FArrayBox& fx,
                                  amrex::FArrayBox& fy,
                                  amrex::FArrayBox& /*fz*/),
                     amrex::GpuArray<amrex::Real,BL_SPACEDIM> dx,
                     const amrex::Real dt)
{
    using namespace amrex;

    static bool inited = false;
    static Real gamma = 1.4;
    static int force_1d = 0;
    static int first_order = 0;
    static int host_update = 1;
    static int debug_state = 0;
    static int test_id = 1;
    static Real x0_prob = 0.5;
    static Real prob_lo_x = 0.0;
    static int ncell_x = -1;
    if (!inited) {
        ParmParse pp("prob");
        pp.query("gamma", gamma);
        pp.query("test", test_id);
        pp.query("x0", x0_prob);
        ParmParse ppadv("adv");
        ppadv.query("force_1d", force_1d);
        ppadv.query("first_order", first_order);
        ppadv.query("host_update", host_update);
        ppadv.query("debug_state", debug_state);

        {
            ParmParse ppgeom("geometry");
            Vector<Real> plo(AMREX_SPACEDIM, 0.0_rt);
            if (ppgeom.countval("prob_lo") >= AMREX_SPACEDIM) {
                ppgeom.queryarr("prob_lo", plo, 0, AMREX_SPACEDIM);
            }
            prob_lo_x = plo[0];
        }
        {
            ParmParse ppamr("amr");
            Vector<int> nc(AMREX_SPACEDIM, -1);
            if (ppamr.countval("n_cell") >= AMREX_SPACEDIM) {
                ppamr.queryarr("n_cell", nc, 0, AMREX_SPACEDIM);
            }
            ncell_x = nc[0];
        }
        inited = true;
    }

    const Real dtdx = dt/dx[0];
    const Real dtdy = (force_1d != 0) ? 0.0_rt : dt/dx[1];

    if (debug_state == 2) {
        // Debug mode: copy input to output directly (should preserve constants)
        auto const Uin_dbg  = statein.const_array();
        auto const Uout_dbg = stateout.array();
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int n=0; n<4; ++n) {
                Uout_dbg(i,j,k,n) = Uin_dbg(i,j,k,n);
            }
        });
        return;
    }
    // Need 2 ghost cells for slope calculation (q(i+/-1)) while slopes live on grow(bx,1)
    const Box gbx = grow(bx,2);
    const Box sbx = grow(bx,1);
    const Box xface = surroundingNodes(bx,0);
    const Box yface = surroundingNodes(bx,1);

    FArrayBox qfab(gbx, NQ);
    FArrayBox sxfab(sbx, NQ);
    FArrayBox syfab(sbx, NQ);
    FArrayBox utmpfab(gbx, NUM_STATE);

    auto const Uin  = statein.const_array();
    auto const Uout = stateout.array();

    auto const q   = qfab.array();
    auto const sx  = sxfab.array();
    auto const sy  = syfab.array();
    auto const Utmp = utmpfab.array();

    auto const Fx = fx.array();
    auto const Fy = fy.array();

    auto compute_fluxes = [&](Array4<Real const> U, Array4<Real> Fxloc, Array4<Real> Fyloc)
    {
        // 1) conserved -> primitive
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real Uc[4] = {U(i,j,k,RHO), U(i,j,k,MX), U(i,j,k,MY), U(i,j,k,E)};
            Real rho,u,v,p;
            cons_to_prim(Uc, gamma, rho, u, v, p);
            q(i,j,k,0) = rho;
            q(i,j,k,1) = u;
            q(i,j,k,2) = v;
            q(i,j,k,3) = p;
        });

        // 2) slopes (minmod) in x and y, or zero for first-order
        ParallelFor(sbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (first_order != 0) {
                for (int n=0; n<NQ; ++n) {
                    sx(i,j,k,n) = 0.0_rt;
                    sy(i,j,k,n) = 0.0_rt;
                }
                return;
            }

            Real sx_loc[NQ];
            Real sy_loc[NQ];
            for (int n=0; n<NQ; ++n) {
                const Real dlx = q(i,j,k,n) - q(i-1,j,k,n);
                const Real drx = q(i+1,j,k,n) - q(i,j,k,n);
                sx_loc[n] = minmod(dlx, drx);

                const Real dly = q(i,j,k,n) - q(i,j-1,k,n);
                const Real dry = q(i,j+1,k,n) - q(i,j,k,n);
                sy_loc[n] = minmod(dly, dry);
            }

            auto limiter = [] AMREX_GPU_DEVICE (Real q0, Real dq, Real qmin) noexcept
            {
                if (dq == 0.0_rt) return 1.0_rt;
                const Real qlo = q0 - 0.5_rt*dq;
                const Real qhi = q0 + 0.5_rt*dq;
                const Real minq = amrex::min(qlo, qhi);
                if (minq < qmin) {
                    return (q0 - qmin) / (q0 - minq + 1.0e-50_rt);
                }
                return 1.0_rt;
            };

            const Real rho0 = q(i,j,k,0);
            const Real p0   = q(i,j,k,3);

            const Real theta_x = amrex::min(limiter(rho0, sx_loc[0], RHO_FLOOR),
                                            limiter(p0,   sx_loc[3], P_FLOOR));
            const Real theta_y = amrex::min(limiter(rho0, sy_loc[0], RHO_FLOOR),
                                            limiter(p0,   sy_loc[3], P_FLOOR));

            for (int n=0; n<NQ; ++n) {
                sx(i,j,k,n) = theta_x * sx_loc[n];
                sy(i,j,k,n) = theta_y * sy_loc[n];
            }
        });

        // 3) HLL fluxes on x-faces
        ParallelFor(xface, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // left cell = i-1, right cell = i
            Real rhoL = q(i-1,j,k,0) + 0.5*sx(i-1,j,k,0);
            Real uL   = q(i-1,j,k,1) + 0.5*sx(i-1,j,k,1);
            Real vL   = q(i-1,j,k,2) + 0.5*sx(i-1,j,k,2);
            Real pL   = q(i-1,j,k,3) + 0.5*sx(i-1,j,k,3);

            Real rhoR = q(i  ,j,k,0) - 0.5*sx(i  ,j,k,0);
            Real uR   = q(i  ,j,k,1) - 0.5*sx(i  ,j,k,1);
            Real vR   = q(i  ,j,k,2) - 0.5*sx(i  ,j,k,2);
            Real pR   = q(i  ,j,k,3) - 0.5*sx(i  ,j,k,3);

            Real F[4];
            hll_flux_x(rhoL,uL,vL,pL, rhoR,uR,vR,pR, gamma, F);
            Fxloc(i,j,k,RHO) = F[0];
            Fxloc(i,j,k,MX)  = F[1];
            Fxloc(i,j,k,MY)  = F[2];
            Fxloc(i,j,k,E)   = F[3];
        });

        // 4) HLL fluxes on y-faces (skip if forcing 1D)
        if (force_1d != 0) {
            ParallelFor(yface, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Fyloc(i,j,k,RHO) = 0.0_rt;
                Fyloc(i,j,k,MX)  = 0.0_rt;
                Fyloc(i,j,k,MY)  = 0.0_rt;
                Fyloc(i,j,k,E)   = 0.0_rt;
            });
        } else {
            ParallelFor(yface, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                // left cell = j-1, right cell = j
                Real rhoL = q(i,j-1,k,0) + 0.5*sy(i,j-1,k,0);
                Real uL   = q(i,j-1,k,1) + 0.5*sy(i,j-1,k,1);
                Real vL   = q(i,j-1,k,2) + 0.5*sy(i,j-1,k,2);
                Real pL   = q(i,j-1,k,3) + 0.5*sy(i,j-1,k,3);

                Real rhoR = q(i,j  ,k,0) - 0.5*sy(i,j  ,k,0);
                Real uR   = q(i,j  ,k,1) - 0.5*sy(i,j  ,k,1);
                Real vR   = q(i,j  ,k,2) - 0.5*sy(i,j  ,k,2);
                Real pR   = q(i,j  ,k,3) - 0.5*sy(i,j  ,k,3);

                Real G[4];
                hll_flux_y(rhoL,uL,vL,pL, rhoR,uR,vR,pR, gamma, G);
                Fyloc(i,j,k,RHO) = G[0];
                Fyloc(i,j,k,MX)  = G[1];
                Fyloc(i,j,k,MY)  = G[2];
                Fyloc(i,j,k,E)   = G[3];
            });
        }
    };

    // Stage 1 fluxes from U^n
    compute_fluxes(Uin, Fx, Fy);
    amrex::Gpu::synchronize();



    if (host_update != 0) {
        // Use the host update path selected by adv.host_update.
        Long pressure_floor_stage1_count = 0;
        Long pressure_floor_stage2_count = 0;
        auto const Uin_h  = statein.const_array();
        auto const Utmp_h = utmpfab.array();
        auto const Uout_h = stateout.array();
        auto const Fx_h = fx.const_array();
        auto const Fy_h = fy.const_array();

        const auto lo = bx.smallEnd();
        const auto hi = bx.bigEnd();
#if (AMREX_SPACEDIM == 3)
        const int klo = lo[2];
        const int khi = hi[2];
#else
        const int klo = 0;
        const int khi = 0;
#endif
        for (int k = klo; k <= khi; ++k) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    const Real fxp0 = Fx_h(i+1,j,k,RHO);
                    const Real fxm0 = Fx_h(i  ,j,k,RHO);
                    const Real fxp1 = Fx_h(i+1,j,k,MX);
                    const Real fxm1 = Fx_h(i  ,j,k,MX);
                    const Real fxp2 = Fx_h(i+1,j,k,MY);
                    const Real fxm2 = Fx_h(i  ,j,k,MY);
                    const Real fxp3 = Fx_h(i+1,j,k,E);
                    const Real fxm3 = Fx_h(i  ,j,k,E);

                    if (force_1d != 0) {
                        Utmp_h(i,j,k,RHO) = Uin_h(i,j,k,RHO) - dtdx*(fxp0-fxm0);
                        Utmp_h(i,j,k,MX)  = Uin_h(i,j,k,MX)  - dtdx*(fxp1-fxm1);
                        Utmp_h(i,j,k,MY)  = Uin_h(i,j,k,MY)  - dtdx*(fxp2-fxm2);
                        Utmp_h(i,j,k,E)   = Uin_h(i,j,k,E)   - dtdx*(fxp3-fxm3);
                    } else {
                        const Real fyp0 = Fy_h(i,j+1,k,RHO);
                        const Real fym0 = Fy_h(i,j  ,k,RHO);
                        const Real fyp1 = Fy_h(i,j+1,k,MX);
                        const Real fym1 = Fy_h(i,j  ,k,MX);
                        const Real fyp2 = Fy_h(i,j+1,k,MY);
                        const Real fym2 = Fy_h(i,j  ,k,MY);
                        const Real fyp3 = Fy_h(i,j+1,k,E);
                        const Real fym3 = Fy_h(i,j  ,k,E);
                        Utmp_h(i,j,k,RHO) = Uin_h(i,j,k,RHO) - dtdx*(fxp0-fxm0) - dtdy*(fyp0-fym0);
                        Utmp_h(i,j,k,MX)  = Uin_h(i,j,k,MX)  - dtdx*(fxp1-fxm1) - dtdy*(fyp1-fym1);
                        Utmp_h(i,j,k,MY)  = Uin_h(i,j,k,MY)  - dtdx*(fxp2-fxm2) - dtdy*(fyp2-fym2);
                        Utmp_h(i,j,k,E)   = Uin_h(i,j,k,E)   - dtdx*(fxp3-fxm3) - dtdy*(fyp3-fym3);
                    }

                    // Positivity fix on U*
                    Real rho = amrex::max(Utmp_h(i,j,k,RHO), RHO_FLOOR);
                    Real mx  = Utmp_h(i,j,k,MX);
                    Real my  = Utmp_h(i,j,k,MY);
                    Real Etot = Utmp_h(i,j,k,E);
                    const Real kinetic = 0.5_rt*(mx*mx + my*my)/rho;
                    Real p = (gamma-1.0_rt)*(Etot - kinetic);
                    if (p < P_FLOOR) {
                        if (AmrLevelAdv::diag_pressure_floor != 0) {
                            ++pressure_floor_stage1_count;
                        }
                        Etot = P_FLOOR/(gamma-1.0_rt) + kinetic;
                    }
                    Utmp_h(i,j,k,RHO) = rho;
                    Utmp_h(i,j,k,E)   = Etot;
                }
            }
        }

        // Stage 2 fluxes from U*
        compute_fluxes(Utmp, Fx, Fy);
        amrex::Gpu::synchronize();

        for (int k = klo; k <= khi; ++k) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    const Real fxp0 = Fx_h(i+1,j,k,RHO);
                    const Real fxm0 = Fx_h(i  ,j,k,RHO);
                    const Real fxp1 = Fx_h(i+1,j,k,MX);
                    const Real fxm1 = Fx_h(i  ,j,k,MX);
                    const Real fxp2 = Fx_h(i+1,j,k,MY);
                    const Real fxm2 = Fx_h(i  ,j,k,MY);
                    const Real fxp3 = Fx_h(i+1,j,k,E);
                    const Real fxm3 = Fx_h(i  ,j,k,E);

                    if (force_1d != 0) {
                        const Real Ustar0 = Utmp_h(i,j,k,RHO) - dtdx*(fxp0-fxm0);
                        const Real Ustar1 = Utmp_h(i,j,k,MX)  - dtdx*(fxp1-fxm1);
                        const Real Ustar2 = Utmp_h(i,j,k,MY)  - dtdx*(fxp2-fxm2);
                        const Real Ustar3 = Utmp_h(i,j,k,E)   - dtdx*(fxp3-fxm3);
                        Uout_h(i,j,k,RHO) = 0.5*(Uin_h(i,j,k,RHO) + Ustar0);
                        Uout_h(i,j,k,MX)  = 0.5*(Uin_h(i,j,k,MX)  + Ustar1);
                        Uout_h(i,j,k,MY)  = 0.5*(Uin_h(i,j,k,MY)  + Ustar2);
                        Uout_h(i,j,k,E)   = 0.5*(Uin_h(i,j,k,E)   + Ustar3);
                    } else {
                        const Real fyp0 = Fy_h(i,j+1,k,RHO);
                        const Real fym0 = Fy_h(i,j  ,k,RHO);
                        const Real fyp1 = Fy_h(i,j+1,k,MX);
                        const Real fym1 = Fy_h(i,j  ,k,MX);
                        const Real fyp2 = Fy_h(i,j+1,k,MY);
                        const Real fym2 = Fy_h(i,j  ,k,MY);
                        const Real fyp3 = Fy_h(i,j+1,k,E);
                        const Real fym3 = Fy_h(i,j  ,k,E);
                        const Real Ustar0 = Utmp_h(i,j,k,RHO) - dtdx*(fxp0-fxm0) - dtdy*(fyp0-fym0);
                        const Real Ustar1 = Utmp_h(i,j,k,MX)  - dtdx*(fxp1-fxm1) - dtdy*(fyp1-fym1);
                        const Real Ustar2 = Utmp_h(i,j,k,MY)  - dtdx*(fxp2-fxm2) - dtdy*(fyp2-fym2);
                        const Real Ustar3 = Utmp_h(i,j,k,E)   - dtdx*(fxp3-fxm3) - dtdy*(fyp3-fym3);
                        Uout_h(i,j,k,RHO) = 0.5*(Uin_h(i,j,k,RHO) + Ustar0);
                        Uout_h(i,j,k,MX)  = 0.5*(Uin_h(i,j,k,MX)  + Ustar1);
                        Uout_h(i,j,k,MY)  = 0.5*(Uin_h(i,j,k,MY)  + Ustar2);
                        Uout_h(i,j,k,E)   = 0.5*(Uin_h(i,j,k,E)   + Ustar3);
                    }

                    // Positivity fix on U^{n+1}
                    Real rho = amrex::max(Uout_h(i,j,k,RHO), RHO_FLOOR);
                    Real mx  = Uout_h(i,j,k,MX);
                    Real my  = Uout_h(i,j,k,MY);
                    Real Etot = Uout_h(i,j,k,E);
                    const Real kinetic = 0.5_rt*(mx*mx + my*my)/rho;
                    Real p = (gamma-1.0_rt)*(Etot - kinetic);
                    if (p < P_FLOOR) {
                        if (AmrLevelAdv::diag_pressure_floor != 0) {
                            ++pressure_floor_stage2_count;
                        }
                        Etot = P_FLOOR/(gamma-1.0_rt) + kinetic;
                    }
                    Uout_h(i,j,k,RHO) = rho;
                    Uout_h(i,j,k,E)   = Etot;
                }
            }
        }

        if (AmrLevelAdv::diag_pressure_floor != 0) {
            AmrLevelAdv::addPressureFloorStage1Count(pressure_floor_stage1_count);
            AmrLevelAdv::addPressureFloorStage2Count(pressure_floor_stage2_count);
        }
    } else {

    // Initialize U* with U^n everywhere (including ghost region),
    // then update only the valid box. This provides ghost values for stage-2 fluxes.
    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        for (int n=0; n<4; ++n) {
            Utmp(i,j,k,n) = Uin(i,j,k,n);
        }
    });

    // U* = U^n - dt * div(F)
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const Real fxp0 = Fx(i+1,j,k,RHO);
        const Real fxm0 = Fx(i  ,j,k,RHO);
        const Real fxp1 = Fx(i+1,j,k,MX);
        const Real fxm1 = Fx(i  ,j,k,MX);
        const Real fxp2 = Fx(i+1,j,k,MY);
        const Real fxm2 = Fx(i  ,j,k,MY);
        const Real fxp3 = Fx(i+1,j,k,E);
        const Real fxm3 = Fx(i  ,j,k,E);

        if (force_1d != 0) {
            Utmp(i,j,k,RHO) = Uin(i,j,k,RHO) - dtdx*(fxp0-fxm0);
            Utmp(i,j,k,MX)  = Uin(i,j,k,MX)  - dtdx*(fxp1-fxm1);
            Utmp(i,j,k,MY)  = Uin(i,j,k,MY)  - dtdx*(fxp2-fxm2);
            Utmp(i,j,k,E)   = Uin(i,j,k,E)   - dtdx*(fxp3-fxm3);
        } else {
            const Real fyp0 = Fy(i,j+1,k,RHO);
            const Real fym0 = Fy(i,j  ,k,RHO);
            const Real fyp1 = Fy(i,j+1,k,MX);
            const Real fym1 = Fy(i,j  ,k,MX);
            const Real fyp2 = Fy(i,j+1,k,MY);
            const Real fym2 = Fy(i,j  ,k,MY);
            const Real fyp3 = Fy(i,j+1,k,E);
            const Real fym3 = Fy(i,j  ,k,E);
            Utmp(i,j,k,RHO) = Uin(i,j,k,RHO) - dtdx*(fxp0-fxm0) - dtdy*(fyp0-fym0);
            Utmp(i,j,k,MX)  = Uin(i,j,k,MX)  - dtdx*(fxp1-fxm1) - dtdy*(fyp1-fym1);
            Utmp(i,j,k,MY)  = Uin(i,j,k,MY)  - dtdx*(fxp2-fxm2) - dtdy*(fyp2-fym2);
            Utmp(i,j,k,E)   = Uin(i,j,k,E)   - dtdx*(fxp3-fxm3) - dtdy*(fyp3-fym3);
        }
    });


    if (AmrLevelAdv::diag_pressure_floor != 0) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
            const Real rho = amrex::max(Utmp(i,j,k,RHO), RHO_FLOOR);
            const Real mx  = Utmp(i,j,k,MX);
            const Real my  = Utmp(i,j,k,MY);
            const Real Ecell = Utmp(i,j,k,E);
            const Real kinetic = 0.5_rt*(mx*mx + my*my)/rho;
            const Real p = (gamma-1.0_rt)*(Ecell - kinetic);
            return { static_cast<Long>(p < P_FLOOR) };
        });
        AmrLevelAdv::addPressureFloorStage1Count(amrex::get<0>(reduce_data.value()));
    }

    // Positivity fix on U*
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real rho = amrex::max(Utmp(i,j,k,RHO), RHO_FLOOR);
        Real mx  = Utmp(i,j,k,MX);
        Real my  = Utmp(i,j,k,MY);
        Real Ecell = Utmp(i,j,k,E);
        const Real kinetic = 0.5_rt*(mx*mx + my*my)/rho;
        Real p = (gamma-1.0_rt)*(Ecell - kinetic);
        if (p < P_FLOOR) {
            Ecell = P_FLOOR/(gamma-1.0_rt) + kinetic;
        }
        Utmp(i,j,k,RHO) = rho;
        Utmp(i,j,k,E)   = Ecell;
    });

    // Fill U* ghost cells by clamping to the valid box before stage-2 fluxes.
    const int ilo = bx.smallEnd(0);
    const int ihi = bx.bigEnd(0);
    const int jlo = bx.smallEnd(1);
    const int jhi = bx.bigEnd(1);
#if (AMREX_SPACEDIM == 3)
    const int klo = bx.smallEnd(2);
    const int khi = bx.bigEnd(2);
#endif
    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
#if (AMREX_SPACEDIM == 3)
        if (i < ilo || i > ihi || j < jlo || j > jhi || k < klo || k > khi) {
            const int ii = amrex::max(ilo, amrex::min(i, ihi));
            const int jj = amrex::max(jlo, amrex::min(j, jhi));
            const int kk = amrex::max(klo, amrex::min(k, khi));
#else
        if (i < ilo || i > ihi || j < jlo || j > jhi) {
            const int ii = amrex::max(ilo, amrex::min(i, ihi));
            const int jj = amrex::max(jlo, amrex::min(j, jhi));
            const int kk = 0;
#endif
            for (int n = 0; n < 4; ++n) {
                Utmp(i,j,k,n) = Utmp(ii,jj,kk,n);
            }
        }
    });


    // Stage 2 fluxes from U*
    compute_fluxes(Utmp, Fx, Fy);


    // U^{n+1} = 0.5*(U^n + U* - dt*div(F(U*)))
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const Real fxp0 = Fx(i+1,j,k,RHO);
        const Real fxm0 = Fx(i  ,j,k,RHO);
        const Real fxp1 = Fx(i+1,j,k,MX);
        const Real fxm1 = Fx(i  ,j,k,MX);
        const Real fxp2 = Fx(i+1,j,k,MY);
        const Real fxm2 = Fx(i  ,j,k,MY);
        const Real fxp3 = Fx(i+1,j,k,E);
        const Real fxm3 = Fx(i  ,j,k,E);

        if (force_1d != 0) {
            const Real Ustar0 = Utmp(i,j,k,RHO) - dtdx*(fxp0-fxm0);
            const Real Ustar1 = Utmp(i,j,k,MX)  - dtdx*(fxp1-fxm1);
            const Real Ustar2 = Utmp(i,j,k,MY)  - dtdx*(fxp2-fxm2);
            const Real Ustar3 = Utmp(i,j,k,E)   - dtdx*(fxp3-fxm3);
            Uout(i,j,k,RHO) = 0.5*(Uin(i,j,k,RHO) + Ustar0);
            Uout(i,j,k,MX)  = 0.5*(Uin(i,j,k,MX)  + Ustar1);
            Uout(i,j,k,MY)  = 0.5*(Uin(i,j,k,MY)  + Ustar2);
            Uout(i,j,k,E)   = 0.5*(Uin(i,j,k,E)   + Ustar3);
        } else {
            const Real fyp0 = Fy(i,j+1,k,RHO);
            const Real fym0 = Fy(i,j  ,k,RHO);
            const Real fyp1 = Fy(i,j+1,k,MX);
            const Real fym1 = Fy(i,j  ,k,MX);
            const Real fyp2 = Fy(i,j+1,k,MY);
            const Real fym2 = Fy(i,j  ,k,MY);
            const Real fyp3 = Fy(i,j+1,k,E);
            const Real fym3 = Fy(i,j  ,k,E);
            const Real Ustar0 = Utmp(i,j,k,RHO) - dtdx*(fxp0-fxm0) - dtdy*(fyp0-fym0);
            const Real Ustar1 = Utmp(i,j,k,MX)  - dtdx*(fxp1-fxm1) - dtdy*(fyp1-fym1);
            const Real Ustar2 = Utmp(i,j,k,MY)  - dtdx*(fxp2-fxm2) - dtdy*(fyp2-fym2);
            const Real Ustar3 = Utmp(i,j,k,E)   - dtdx*(fxp3-fxm3) - dtdy*(fyp3-fym3);
            Uout(i,j,k,RHO) = 0.5*(Uin(i,j,k,RHO) + Ustar0);
            Uout(i,j,k,MX)  = 0.5*(Uin(i,j,k,MX)  + Ustar1);
            Uout(i,j,k,MY)  = 0.5*(Uin(i,j,k,MY)  + Ustar2);
            Uout(i,j,k,E)   = 0.5*(Uin(i,j,k,E)   + Ustar3);
        }
    });
    amrex::Gpu::synchronize();



    if (AmrLevelAdv::diag_pressure_floor != 0) {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
            const Real rho = amrex::max(Uout(i,j,k,RHO), RHO_FLOOR);
            const Real mx  = Uout(i,j,k,MX);
            const Real my  = Uout(i,j,k,MY);
            const Real Ecell = Uout(i,j,k,E);
            const Real kinetic = 0.5_rt*(mx*mx + my*my)/rho;
            const Real p = (gamma-1.0_rt)*(Ecell - kinetic);
            return { static_cast<Long>(p < P_FLOOR) };
        });
        AmrLevelAdv::addPressureFloorStage2Count(amrex::get<0>(reduce_data.value()));
    }

    // Positivity fix on U^{n+1}
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real rho = amrex::max(Uout(i,j,k,RHO), RHO_FLOOR);
        const Real mx  = Uout(i,j,k,MX);
        const Real my  = Uout(i,j,k,MY);
        Real Ecell = Uout(i,j,k,E);
        const Real kinetic = 0.5_rt*(mx*mx + my*my)/rho;
        Real p = (gamma-1.0_rt)*(Ecell - kinetic);
        if (p < P_FLOOR) {
            Ecell = P_FLOOR/(gamma-1.0_rt) + kinetic;
        }
        Uout(i,j,k,RHO) = rho;
        Uout(i,j,k,E)   = Ecell;
    });
    amrex::Gpu::synchronize();



    }

    // Copy fluxes into Flux MultiFab (already in Fx/Fy arrays)
    // Scale by face area for reflux
    ParallelFor(
        AMREX_D_DECL(xface, yface, Box()),
        AMREX_D_DECL([=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                     {
                         Fx(i,j,k,RHO) *= dt*dx[1];
                         Fx(i,j,k,MX)  *= dt*dx[1];
                         Fx(i,j,k,MY)  *= dt*dx[1];
                         Fx(i,j,k,E)   *= dt*dx[1];
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                     {
                         Fy(i,j,k,RHO) *= dt*dx[0];
                         Fy(i,j,k,MX)  *= dt*dx[0];
                         Fy(i,j,k,MY)  *= dt*dx[0];
                         Fy(i,j,k,E)   *= dt*dx[0];
                     },
                     [=] AMREX_GPU_DEVICE (int, int, int) noexcept {}));

    amrex::ignore_unused(nbx);
}
