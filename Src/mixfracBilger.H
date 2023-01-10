#ifndef MIXFRACBILGER_H
#define MIXFRACBILGER_H

#include <mechanism.H>

struct mixFracBilgerData {
   amrex::Array<amrex::Real, 4> Beta_mix;
   amrex::Array<amrex::Real, NUM_SPECIES> spec_Bilger_fact;
   amrex::Real Zfu = -1.0;
   amrex::Real Zox = -1.0;
};
#endif
