// Modified for the Euler AMR project - the portfolio version
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Reduce.H>
#include <AMReX_Math.H>
#include <AMReX_BCUtil.H>

#include "AmrLevelAdv.H"
#include "Prob.H"
#include "Kernels.H"

using namespace amrex;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9;
int      AmrLevelAdv::do_reflux       = 1;
int      AmrLevelAdv::diag_timers     = 0;
int      AmrLevelAdv::diag_pressure_floor = 0;
Real     AmrLevelAdv::diag_fillpatch_time = 0.0_rt;
Real     AmrLevelAdv::diag_advance_time = 0.0_rt;
Real     AmrLevelAdv::diag_reflux_time = 0.0_rt;
Real     AmrLevelAdv::diag_tagging_time = 0.0_rt;
Long     AmrLevelAdv::diag_pressure_floor_stage1 = 0;
Long     AmrLevelAdv::diag_pressure_floor_stage2 = 0;

// Number of conserved variables for Euler (rho, momx, momy, [momz], E)
#if (AMREX_SPACEDIM == 2)
int      AmrLevelAdv::NUM_STATE       = 4;
#elif (AMREX_SPACEDIM == 3)
int      AmrLevelAdv::NUM_STATE       = 5;
#else
int      AmrLevelAdv::NUM_STATE       = 3; // fallback if 1D is enabled elsewhere
#endif
int      AmrLevelAdv::NUM_GROW        = 3;  // number of ghost cells

ProbParm* AmrLevelAdv::h_prob_parm = nullptr;
ProbParm* AmrLevelAdv::d_prob_parm = nullptr;

int      AmrLevelAdv::max_phierr_lev  = -1;
int      AmrLevelAdv::max_phigrad_lev = -1;

Vector<Real> AmrLevelAdv::phierr;
Vector<Real> AmrLevelAdv::phigrad;
Vector<BCRec> AmrLevelAdv::bcs;

#ifdef AMREX_PARTICLES
std::unique_ptr<AmrTracerParticleContainer> AmrLevelAdv::TracerPC =  nullptr;
int AmrLevelAdv::do_tracers                       =  0;
#endif

AmrLevelAdv::AmrLevelAdv ()
{
    flux_reg = nullptr;
}

AmrLevelAdv::AmrLevelAdv (Amr&            papa,
                          int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    flux_reg = nullptr;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

AmrLevelAdv::~AmrLevelAdv ()
{
    delete flux_reg;
}

void
AmrLevelAdv::restart (Amr&          papa,
                      std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

void
AmrLevelAdv::checkPoint (const std::string& dir,
                         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
#ifdef AMREX_PARTICLES
  if (do_tracers && level == 0) {
    TracerPC->Checkpoint(dir, "Tracer", true);
  }
#endif
}

void
AmrLevelAdv::writePlotFile (const std::string& dir,
                             std::ostream&      os,
                            VisMF::How         how)
{

    AmrLevel::writePlotFile (dir,os,how);

#ifdef AMREX_PARTICLES
    if (do_tracers && level == 0) {
      TracerPC->Checkpoint(dir, "Tracer", true);
    }
#endif
}

void
AmrLevelAdv::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Initialize struct containing problem-specific variables
    h_prob_parm = new ProbParm{};
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));

    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
                           &cell_cons_interp);

    Geometry const* gg = AMReX::top()->getDefaultGeometry();
    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        if (gg->isPeriodic(i)) {
            lo_bc[i] = hi_bc[i] = BCType::int_dir;  // periodic
        } else {
            lo_bc[i] = hi_bc[i] = BCType::foextrap; // outflow
        }
    }

    BCRec bc(lo_bc, hi_bc);

    StateDescriptor::BndryFunc bndryfunc(nullfill);
    bndryfunc.setRunOnGPU(true);

    // Register conserved Euler variables
#if (AMREX_SPACEDIM == 2)
    const std::vector<std::string> names = {"rho","momx","momy","E"};
#else
    const std::vector<std::string> names = {"rho","momx","momy","momz","E"};
#endif
    bcs.resize(NUM_STATE);
    for (int n = 0; n < NUM_STATE; ++n) {
        desc_lst.setComponent(Phi_Type, n, names[n], bc, bndryfunc);
        bcs[n] = bc;
    }
}

void
AmrLevelAdv::variableCleanUp ()
{
    desc_lst.clear();
#ifdef AMREX_PARTICLES
    TracerPC.reset();
#endif

    // Delete structs containing problem-specific parameters
    delete h_prob_parm;
    The_Arena()->free(d_prob_parm);
}

void
AmrLevelAdv::resetDiagnostics ()
{
    diag_fillpatch_time = 0.0_rt;
    diag_advance_time = 0.0_rt;
    diag_reflux_time = 0.0_rt;
    diag_tagging_time = 0.0_rt;
    diag_pressure_floor_stage1 = 0;
    diag_pressure_floor_stage2 = 0;
}

void
AmrLevelAdv::addFillPatchTime (Real dt)
{
    diag_fillpatch_time += dt;
}

void
AmrLevelAdv::addAdvanceTime (Real dt)
{
    diag_advance_time += dt;
}

void
AmrLevelAdv::addRefluxTime (Real dt)
{
    diag_reflux_time += dt;
}

void
AmrLevelAdv::addTaggingTime (Real dt)
{
    diag_tagging_time += dt;
}

void
AmrLevelAdv::addPressureFloorStage1Count (Long count)
{
    diag_pressure_floor_stage1 += count;
}

void
AmrLevelAdv::addPressureFloorStage2Count (Long count)
{
    diag_pressure_floor_stage2 += count;
}

void
AmrLevelAdv::reportDiagnostics (Real total_runtime)
{
    if (diag_timers == 0 && diag_pressure_floor == 0) {
        return;
    }

    const int io_proc = ParallelDescriptor::IOProcessorNumber();

    if (diag_timers != 0) {
        Real fillpatch = diag_fillpatch_time;
        Real advance = diag_advance_time;
        Real reflux = diag_reflux_time;
        Real tagging = diag_tagging_time;
        ParallelDescriptor::ReduceRealMax(fillpatch, io_proc);
        ParallelDescriptor::ReduceRealMax(advance, io_proc);
        ParallelDescriptor::ReduceRealMax(reflux, io_proc);
        ParallelDescriptor::ReduceRealMax(tagging, io_proc);

        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print()
                << "Diag timing summary (max over MPI ranks): "
                << "fillpatch=" << fillpatch
                << " advance=" << advance
                << " reflux=" << reflux
                << " tagging=" << tagging
                << " total_runtime=" << total_runtime
                << std::endl;
            if (total_runtime > 0.0_rt) {
                amrex::Print()
                    << "Diag timing fractions (% of total runtime): "
                    << "fillpatch=" << 100.0_rt * fillpatch / total_runtime
                    << " advance=" << 100.0_rt * advance / total_runtime
                    << " reflux=" << 100.0_rt * reflux / total_runtime
                    << " tagging=" << 100.0_rt * tagging / total_runtime
                    << std::endl;
            }
        }
    }

    if (diag_pressure_floor != 0) {
        Long stage1 = diag_pressure_floor_stage1;
        Long stage2 = diag_pressure_floor_stage2;
        ParallelDescriptor::ReduceLongSum(stage1, io_proc);
        ParallelDescriptor::ReduceLongSum(stage2, io_proc);

        if (ParallelDescriptor::IOProcessor()) {
            amrex::Print()
                << "Diag pressure-floor triggers (sum over MPI ranks): "
                << "stage1=" << stage1
                << " stage2=" << stage2
                << " total=" << (stage1 + stage2)
                << std::endl;
        }
    }
}

/**
 * Initialize grid data at problem start-up.
 */
void
AmrLevelAdv::initData ()
{
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

    // C++ initial conditions (problem-specific)
    initdata(S_new, geom, cur_time);

    // Initialize old data to match new data so FillPatch has valid time levels
    state[Phi_Type].allocOldData();
    state[Phi_Type].oldData().ParallelCopy(S_new, 0, 0, NUM_STATE, 0, 0);
    state[Phi_Type].setOldTimeLevel(cur_time);

#ifdef AMREX_PARTICLES
    init_particles();
#endif

    if (verbose) {
        amrex::Print() << "Done initializing the level " << level
                       << " data " << std::endl;
    }
}

/**
 * Initialize data on this level from another AmrLevelAdv (during regrid).
 */
void
AmrLevelAdv::init (AmrLevel &old)
{
    auto* oldlev = static_cast<AmrLevelAdv*>(&old);

    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

/**
 * Initialize data on this level after regridding if old level did not previously exist
 */
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

/**
 * Advance grids at this level in time.
 */
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  /*ncycle*/)
{
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    MultiFab& S_new = get_new_data(Phi_Type);

    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);

    GpuArray<Real,BL_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,BL_SPACEDIM> prob_lo = geom.ProbLoArray();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;

    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
        fine = &getFluxReg(level+1);
        fine->setVal(0.0);
    }

    if (do_reflux && level > 0) {
        current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM];

    if (do_reflux)
    {
        for (int j = 0; j < BL_SPACEDIM; j++)
        {
            BoxArray ba = S_new.boxArray();
            ba.surroundingNodes(j);
            fluxes[j].define(ba, dmap, NUM_STATE, 0);
        }
    }

    // State with ghost cells
    MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
    const auto fillpatch_start = (diag_timers != 0) ? amrex::second() : 0.0_rt;
    FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    Sborder.FillBoundary(geom.periodicity());
    FillDomainBoundary(Sborder, geom, bcs);
    if (diag_timers != 0) {
        amrex::Gpu::synchronize();
        addFillPatchTime(amrex::second() - fillpatch_start);
    }


    // Dummy face velocity MultiFabs (not used by Euler update)
    MultiFab Umac[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      BoxArray ba = S_new.boxArray();
      ba.surroundingNodes(i);
      Umac[i].define(ba, dmap, 1, iteration);
    }

    const auto advance_start = (diag_timers != 0) ? amrex::second() : 0.0_rt;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox fluxfab[AMREX_SPACEDIM];
        FArrayBox* flux[AMREX_SPACEDIM];
        FArrayBox* uface[AMREX_SPACEDIM];

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            // Set up tileboxes and nodal tileboxes
            const Box& bx = mfi.tilebox();
            GpuArray<Box,BL_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);,
                         nbx[1] = mfi.nodaltilebox(1);,
                         nbx[2] = mfi.nodaltilebox(2));

            // Grab fab pointers from state multifabs
            const FArrayBox& statein = Sborder[mfi];
            FArrayBox& stateout      =   S_new[mfi];

            for (int i = 0; i < BL_SPACEDIM ; i++) {
#ifdef AMREX_USE_GPU
                // No tiling on GPU.
                flux[i] = &(fluxes[i][mfi]);
                uface[i] = &(Umac[i][mfi]); // unused in Euler update
#else
                const Box& bxtmp = amrex::surroundingNodes(bx,i);
                fluxfab[i].resize(bxtmp,NUM_STATE);
                flux[i] = &(fluxfab[i]);
                uface[i] = &(Umac[i][mfi]); // unused in Euler update
#endif
            }

            // Advect. See Adv.cpp for implementation.
            advect(time, bx, nbx, statein, stateout,
                   AMREX_D_DECL(*uface[0], *uface[1], *uface[2]),
                   AMREX_D_DECL(*flux[0],  *flux[1],  *flux[2]),
                   dx, dt);

#ifndef AMREX_USE_GPU
            if (do_reflux) {
                for (int i = 0; i < BL_SPACEDIM ; i++)
                    fluxes[i][mfi].copy(*flux[i],mfi.nodaltilebox(i));
            }
#endif
        }
    }
    if (diag_timers != 0) {
        amrex::Gpu::synchronize();
        addAdvanceTime(amrex::second() - advance_start);
    }


    if (do_reflux) {
        if (current) {
            for (int i = 0; i < BL_SPACEDIM ; i++)
                current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
        }
        if (fine) {
            for (int i = 0; i < BL_SPACEDIM ; i++)
                fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
        }
    }

#ifdef AMREX_PARTICLES
    if (TracerPC) {
      TracerPC->AdvectWithUmac(Umac, level, dt);
    }
#endif

    return dt;
}

/**
 * Estimate time step.
 */
Real
AmrLevelAdv::estTimeStep (Real)
{
    Real dt_est  = 1.0e+20;

    static bool inited = false;
    static Real gamma = 1.4;
    if (!inited) {
        ParmParse pp("prob");
        pp.query("gamma", gamma);
        inited = true;
    }

    GpuArray<Real,BL_SPACEDIM> dx = geom.CellSizeArray();
    const MultiFab& S_new = get_new_data(Phi_Type);

    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& U = S_new.const_array(mfi);

        amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
        amrex::ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            const Real rho = U(i,j,k,0);
            const Real u   = U(i,j,k,1)/rho;
            const Real v   = U(i,j,k,2)/rho;
            const Real E   = U(i,j,k,3);
            const Real p   = (gamma-1.0_rt)*(E - 0.5_rt*rho*(u*u+v*v));
            const Real c   = std::sqrt(gamma*p/rho);
            const Real ax  = amrex::Math::abs(u) + c;
            const Real ay  = amrex::Math::abs(v) + c;
            const Real dt_cell = 1.0_rt / (ax/dx[0] + ay/dx[1]);
            return {dt_cell};
        });

        dt_est = std::min(dt_est, amrex::get<0>(reduce_data.value()));
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    if (verbose) {
        amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level
                       << ":  dt_est = " << dt_est << std::endl;
    }

    return dt_est;
}

/**
 * Compute initial time step.
 */
Real
AmrLevelAdv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

/**
 * Compute initial `dt'.
 */
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
                               int                   /*sub_cycle*/,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& /*ref_ratio*/,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

/**
 * Compute new `dt'.
 */
void
AmrLevelAdv::computeNewDt (int                   finest_level,
                           int                   /*sub_cycle*/,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& /*ref_ratio*/,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        AmrLevelAdv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1)
    {
        //
        // Limit dt's by pre-regrid dt
        //
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],dt_level[i]);
        }
    }
    else
    {
        //
        // Limit dt's by change_max * old dt
        //
        static Real change_max = 1.1;
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        }
    }

    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

/**
 * Do work after timestep().
 */
void
AmrLevelAdv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();

#ifdef AMREX_PARTICLES
    if (TracerPC)
      {
        const int ncycle = parent->nCycle(level);

        if (iteration < ncycle || level == 0)
          {
            int ngrow = (level == 0) ? 0 : iteration;

            TracerPC->Redistribute(level, TracerPC->finestLevel(), ngrow);
          }
      }
#endif
}

/**
 * Do work after regrid().
 */
void
AmrLevelAdv::post_regrid (int lbase, int /*new_finest*/) {
#ifdef AMREX_PARTICLES
  if (TracerPC && level == lbase) {
      TracerPC->Redistribute(lbase);
  }
#else
  amrex::ignore_unused(lbase);
#endif
}

/**
 * Do work after a restart().
 */
void
AmrLevelAdv::post_restart()
{
#ifdef AMREX_PARTICLES
    if (do_tracers && level == 0) {
      BL_ASSERT(TracerPC == 0);
      TracerPC = std::make_unique<AmrTracerParticleContainer>(parent);
      TracerPC->Restart(parent->theRestartFile(), "Tracer");
    }
#endif
}

/**
 * Do work after init().
 */
void
AmrLevelAdv::post_init (Real /*stop_time*/)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

/**
 * Error estimation for regridding.
 */
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
                       int          /*clearval*/,
                       int          /*tagval*/,
                       Real         /*time*/,
                       int          /*n_error_buf*/,
                       int          /*ngrow*/)
{
    const auto tagging_start = (diag_timers != 0) ? amrex::second() : 0.0_rt;
    MultiFab& S_new = get_new_data(Phi_Type);

    // Properly fill patches and ghost cells for phi gradient check.
    MultiFab phitmp;
    if (level < max_phigrad_lev) {
        const Real cur_time = state[Phi_Type].curTime();
        phitmp.define(S_new.boxArray(), S_new.DistributionMap(), NUM_STATE, 1);
        FillPatch(*this, phitmp, 1, cur_time, Phi_Type, 0, NUM_STATE);
        phitmp.FillBoundary(geom.periodicity());
        FillDomainBoundary(phitmp, geom, bcs);
    }
    MultiFab const& phi = (level < max_phigrad_lev) ? phitmp : S_new;

    const char   tagval = TagBox::SET;
    // const char clearval = TagBox::CLEAR;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& tilebx  = mfi.tilebox();
            const auto phiarr  = phi.array(mfi);
            auto       tagarr  = tags.array(mfi);

            // Tag cells with high phi.
            if (level < max_phierr_lev) {
                const Real phierr_lev  = phierr[level];
                amrex::ParallelFor(tilebx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    state_error(i, j, k, tagarr, phiarr, phierr_lev, tagval);
                });
            }

            // Tag cells with high phi gradient.
            if (level < max_phigrad_lev) {
                const Real phigrad_lev = phigrad[level];
                amrex::ParallelFor(tilebx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    grad_error(i, j, k, tagarr, phiarr, phigrad_lev, tagval);
                });
            }
        }
    }

    if (diag_timers != 0) {
        amrex::Gpu::synchronize();
        addTaggingTime(amrex::second() - tagging_start);
    }
}

/**
 * Read parameters from input file.
 */
void
AmrLevelAdv::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);
    pp.query("diag_timers",diag_timers);
    pp.query("diag_pressure_floor",diag_pressure_floor);

    Geometry const* gg = AMReX::top()->getDefaultGeometry();

    // The solver assumes Cartesian geometry.
    if (! gg->IsCartesian()) {
        amrex::Abort("Please set geom.coord_sys = 0");
    }

    // Allow non-periodic boundaries for shock-tube tests.

#ifdef AMREX_PARTICLES
    pp.query("do_tracers", do_tracers);
#endif

    // Read tagging parameters from tagging block in the input file.
    // See Src_nd/Tagging_params.cpp for the function implementation.
    get_tagging_params();
}

void
AmrLevelAdv::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const auto strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        auto      end    = amrex::second() - strt;

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        amrex::Print() << "AmrLevelAdv::reflux() at level " << level
                       << " : time = " << end << std::endl;
    }

    if (diag_timers != 0) {
        addRefluxTime(amrex::second() - strt);
    }
}

void
AmrLevelAdv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelAdv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);

    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

#ifdef AMREX_PARTICLES
void
AmrLevelAdv::init_particles ()
{
  if (do_tracers && level == 0)
    {
      BL_ASSERT(TracerPC == nullptr);

      TracerPC = std::make_unique<AmrTracerParticleContainer>(parent);

      AmrTracerParticleContainer::ParticleInitData pdata = {{AMREX_D_DECL(0.0, 0.0, 0.0)},{},{},{}};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
    }
}
#endif
