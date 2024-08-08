/************************************************************************

  smfcat

  This is a particularly simple filter, much like the Unix "cat" filter.
  It reads all the specified inputs, concatenates them into a single
  model, and outputs the results in a canonical form.

  Copyright (C) 1998 Michael Garland, All Rights Reserved.
  
  $Id: smfcat.cxx,v 1.1.1.1 2006/09/20 01:42:05 marc Exp $

 ************************************************************************/

#include <stdmix.h>
#include "cmdline.h"

main(int argc, char *argv[])
{
    MxStdModel *m = process_cmdline(argc, argv);
    if( !m ) return 0;

    output_final_model(m);
    delete m;
    return 0;
}
