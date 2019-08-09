#include <util.H>
#include <util_F.H>

#include <AMReX_Print.H>

using namespace amrex;

static std::string
decodeStringFromFortran(const Vector<int>& coded,
                        int                length)
{
    std::string result;
    for (int i = 0; i < length; ++i)
        result += coded[i];
    return result;
}

namespace analysis_util {

  int
  NumSpecies()
  {
    return num_species();
  }
  
  int
  NumElements()
  {
    return num_elements();
  }
  
  Vector<std::string>
  GetSpecNames ()
  {
    int nspecies = analysis_util::NumSpecies();
    Vector<std::string> spec_names(nspecies);
    int Lmax = max_spec_namelen();
    Vector<int> coded_name(Lmax);
    for (int i = 0; i < nspecies; i++)
    {
      int len = Lmax;
      get_spec_names(coded_name.dataPtr(),&i,&len);
      spec_names[i] = decodeStringFromFortran(coded_name,len);
    }
    return spec_names;
  }

  Vector<std::string>
  GetElemNames ()
  {
    int nelements = analysis_util::NumElements();
    Vector<std::string> elem_names(nelements);
    int Lmax = max_elem_namelen();
    Vector<int> coded_name(Lmax);
    for (int i = 0; i < nelements; i++)
    {
      int len = Lmax;
      get_elem_names(coded_name.dataPtr(),&i,&len);
      elem_names[i] = decodeStringFromFortran(coded_name,len);
    }
    return elem_names;
  }

}


