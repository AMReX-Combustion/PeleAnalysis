#ifndef MXGLUTILS_INCLUDED // -*- C++ -*-
#define MXGLUTILS_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  MxGLUtils

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxGLUtils.h,v 1.1.1.1 2006/09/20 01:42:05 marc Exp $

 ************************************************************************/

#include "MxGL.h"

extern GLenum mxgl_matrix_for_mode(GLenum mode);

extern int mxgl_unproject_pixel(int *pixel, double *world, double z=0.0);

extern void mxgl_begin_pick(int *where, double rad, GLuint *buffer, int size);
extern unsigned int mxgl_complete_pick(GLuint *buffer);

// MXGLUTILS_INCLUDED
#endif
