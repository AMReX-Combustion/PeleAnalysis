#ifndef MXSTACK_INCLUDED // -*- C++ -*-
#define MXSTACK_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  This provides a very simple typed-access stack class.  It's really
  just a convenience wrapper for the underlying MxDynBlock class.

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxStack.h,v 1.1.1.1 2006/09/20 01:42:05 marc Exp $

 ************************************************************************/

#include "MxDynBlock.h"

template<class T>
class MxStack : private MxDynBlock<T>
{
public:
    MxStack(unsigned int n) : MxDynBlock<T>(n)
	{ }
    MxStack(const T& val, unsigned int n) : MxDynBlock<T>(n)
	{ push(val); }

    T&       top()       { return MxDynBlock<T>::last(); }
    const T& top() const { return MxDynBlock<T>::last(); }

    bool is_empty() { return MxDynBlock<T>::length()==0; }

    T& pop() { return MxDynBlock<T>::drop(); }

    void push(const T& val) { MxDynBlock<T>::add(val); }
    //
    // NOTE: In this code, it is *crucial* that we do the add() and
    //       assignment in separate steps.  The obvious alternative
    //       is something like { add(top()); }.  But this is subtly
    //       broken!  The top() will grab a pointer into the block,
    //       but the add() may reallocate the block before doing the
    //       assignment.  Thus, the pointer will become invalid.
    void push() { MxDynBlock<T>::add(); top() = (*this)[MxDynBlock<T>::length()-2]; }
};

// MXSTACK_INCLUDED
#endif
