#ifndef util_F_H
#define util_F_H
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
  
  void init_mech();

  void get_spec_names(int* spec_names, int* ispec, int* len);

#ifdef __cplusplus
};
#endif

#endif
  
