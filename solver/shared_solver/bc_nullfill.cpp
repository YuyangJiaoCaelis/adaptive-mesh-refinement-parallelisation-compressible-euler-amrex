#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

struct OutflowFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
            // First-order extrapolation (outflow): copy nearest interior cell.
            const Box& dom = geom.Domain();
            const int i = iv[0];
            const int j = iv[1];
            const int k = iv[2];

            const int ii = amrex::min(amrex::max(i, dom.smallEnd(0)), dom.bigEnd(0));
#if (AMREX_SPACEDIM >= 2)
            const int jj = amrex::min(amrex::max(j, dom.smallEnd(1)), dom.bigEnd(1));
#else
            const int jj = j;
#endif
#if (AMREX_SPACEDIM == 3)
            const int kk = amrex::min(amrex::max(k, dom.smallEnd(2)), dom.bigEnd(2));
#else
            const int kk = k;
#endif

            for (int n = 0; n < numcomp; ++n) {
                dest(i,j,k,dcomp+n) = dest(ii,jj,kk,dcomp+n);
            }
        }
};

void nullfill (Box const& bx, FArrayBox& data,
               const int dcomp, const int numcomp,
               Geometry const& geom, const Real time,
               const Vector<BCRec>& bcr, const int bcomp,
               const int scomp)
{
    GpuBndryFuncFab<OutflowFill> gpu_bndry_func(OutflowFill{});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}
