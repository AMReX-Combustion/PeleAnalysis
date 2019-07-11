#include <util.H>
#include <util_F.H>
#include <mechanism.h>

using namespace amrex;

namespace analysis_util {

  Vector<std::string>
  GetSpecNames ()
  {
    Vector<std::string> spec_names(NUM_SPECIES);
    for (int i = 0; i < NUM_SPECIES; i++)
    {
      int len = 20;
      Vector<int> int_spec_names(len);
      // This call return the actual length of each string in "len" 
      get_spec_names(int_spec_names.dataPtr(),&i,&len);
      char char_spec_names[len+1];
      for (int j = 0; j < len; j++) 
        char_spec_names[j] = int_spec_names[j];
      char_spec_names[len] = '\0';
      spec_names[i] = std::string(char_spec_names);
    }
    return spec_names;
  }

}


