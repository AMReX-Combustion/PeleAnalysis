#ifndef QSLIM_INCLUDED // -*- C++ -*-
#define QSLIM_INCLUDED
#if !defined(__GNUC__)
#  pragma once
#endif

/************************************************************************

  QSlim header file

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: qslim.h,v 1.1.1.1 2006/09/20 01:42:04 marc Exp $

 ************************************************************************/

#include <MxQSlim.h>
#include <MxSMF.h>

#define QSLIM_VERSION 2100
#define QSLIM_VERSION_STRING "2.1"

typedef MxDynBlock<MxPairContraction> QSlimLog;

enum OutputFormat { SMF, PM, MMF, LOG, IV, VRML, MEF };
enum InputFormat { iSMF, iMEF };

////////////////////////////////////////////////////////////////////////
//
// Globally visible (configuration) variables
//

extern unsigned int face_target;
extern bool will_use_fslim;
extern int placement_policy;
extern double boundary_weight;
extern int weighting_policy;
extern bool will_record_history;
extern double compactness_ratio;
extern double meshing_penalty;
extern bool will_join_only;
extern bool be_quiet;
extern OutputFormat output_format;
extern InputFormat input_format;
extern char *output_filename;

extern MxSMFReader *smf;
extern MxStdModel *m;
extern MxStdModel *m_orig;
extern MxQSlim *slim;
extern MxEdgeQSlim *eslim;
extern MxFaceQSlim *fslim;
extern QSlimLog *history;
extern MxDynBlock<MxEdge> *target_edges;

extern bool (*unparsed_hook)(char *, int, char*[], MxStdModel&); 

////////////////////////////////////////////////////////////////////////
//
// Command line parsing and application initialization
//

extern void process_cmdline(int argc, char **argv);
extern void slim_init();
extern void slim_cleanup();
extern void input_file(const char *);
extern void defer_file_inclusion(char *);
extern void include_deferred_files();
extern void slim_history_callback(const MxPairContraction&,float);
extern bool select_input_format(const char *);

////////////////////////////////////////////////////////////////////////
//
// Output routines
//
extern bool select_output_format(const char *);
extern void output_preamble();
extern void output_current_model();
extern void output_final_model();

extern void output_iv(ostream&);
extern void output_vrml(ostream&);
extern void output_regressive_mmf(ostream&);
extern void output_regressive_log(ostream&);
extern void output_progressive_pm(ostream&);

////////////////////////////////////////////////////////////////////////
//
// Other relevant things
//
extern const char *slim_copyright_notice;
extern const char *slim_version_string;
extern void slim_print_banner(ostream&);

// QSLIM_INCLUDED
#endif
