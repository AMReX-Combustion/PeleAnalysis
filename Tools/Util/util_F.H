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
  int num_reactions();
  
  int max_spec_namelen();
  int max_elem_namelen();

  void get_spec_names(int* spec_names, int* ispec, int* len);
  void get_elem_names(int* elem_names, int* ielem, int* len);

  void get_spec_nu_for_rn(int* ids, int* nu, int* len,int* ireac);
  void get_rns_with_X_on_L(int* rns, int * len, int * i);
  void get_rns_with_X_on_R(int* rns, int * len, int * i);

  void get_rns_map(int* rmap);
  void get_elt_comp(int* ispec, int* comp);

#ifdef __cplusplus
};
#endif

#endif
  
