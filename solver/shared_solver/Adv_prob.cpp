// Modified for the Euler AMR project - public project version
#include <AMReX_REAL.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex::Real* /*problo*/,
                         const amrex::Real* /*probhi*/)
    {
        // No separate probin initialisation is required.
    }
}
