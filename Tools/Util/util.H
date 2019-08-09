#ifndef util_H
#define util_H

#include <AMReX_Vector.H>

namespace analysis_util
{
  int NumElements();
  int NumSpecies();
  
  amrex::Vector<std::string> GetElemNames();
  amrex::Vector<std::string> GetSpecNames();
}
#endif
