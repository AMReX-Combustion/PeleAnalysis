#ifndef MXCMDPARSER_INCLUDED // -*- C++ -*-
#define MXCMDPARSER_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  MxCmdParser

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: MxCmdParser.h,v 1.1.1.1 2006/09/20 01:42:05 marc Exp $

 ************************************************************************/

#include "MxDynBlock.h"
#include "MxAsp.h"

typedef MxDynBlock<char *> MxCmdPhrase;

class MxCmd
{
public:
    char *op;
    MxDynBlock<MxCmdPhrase> phrases;

    MxCmd(int N) : phrases(N) { op=NULL; }
};

class MxCmdParser
{
private:
    MxCmd cmd;
    MxAspStore store;

public:
    bool will_ignore_unknown;

public:
    MxCmdParser();

    MxAspStore *asp_store() { return &store; }

    void parse_line(char *, void *closure=NULL);

    virtual bool execute_command(const MxCmd& cmd, void *closure=NULL);
};

// MXCMDPARSER_INCLUDED
#endif
