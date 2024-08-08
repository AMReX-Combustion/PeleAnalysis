#ifndef MXMAT3_INCLUDED // -*- C++ -*-
#define MXMAT3_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  3x3 Matrix class

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxMat3.h,v 1.1.1.1 2006/09/20 01:42:05 marc Exp $

 ************************************************************************/

#include <gfx/mat3.h>

extern bool jacobi(const Mat3& m, Vec3& vals, Vec3 vecs[3]);
extern bool jacobi(const Mat3& m, double *vals, double *vecs);

extern bool fast_jacobi(const Mat3& m, Vec3& vals, Vec3 vecs[3]);

// MXMAT3_INCLUDED
#endif
