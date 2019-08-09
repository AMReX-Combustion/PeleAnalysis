#ifndef util_F_H
#define util_F_H
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
  
  void init_mech();

  int num_species();
  int num_elements();
  int max_spec_namelen();
  int max_elem_namelen();

  void get_spec_names(int* spec_names, int* ispec, int* len);
  void get_elem_names(int* elem_names, int* ielem, int* len);

#ifdef __cplusplus
};
#endif

#endif
  
