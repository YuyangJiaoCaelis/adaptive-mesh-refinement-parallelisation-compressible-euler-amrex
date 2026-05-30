#ifndef PROB_H_
#define PROB_H_

#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

// Initialize conserved Euler state U = [rho, momx, momy, (momz), E]
void initdata(amrex::MultiFab& state, amrex::Geometry const& geom, amrex::Real time);

#endif
