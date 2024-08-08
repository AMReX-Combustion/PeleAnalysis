#ifndef GFXMFC_INCLUDED // -*- C++ -*-
#define GFXMFC_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  Support code for using MFC.  At the moment, this just makes sure that
  we include the right headers.

  $Id: mfc.h,v 1.1.1.1 2006/09/20 01:42:04 marc Exp $

 ************************************************************************/

#ifndef _MBCS
#define _MBCS
#endif

#ifndef _AFXDLL
#define _AFXDLL
#endif

#include <afxwin.h>
#include <afxext.h>

#include "wintools.h"

// GFXMFC_INCLUDED
#endif
