/************************************************************************

  Main QSlim file.

  Copyright (C) 1998 Michael Garland.  See "COPYING.txt" for details.
  
  $Id: main.cpp,v 1.2 2006/09/20 17:41:08 marc Exp $

 ************************************************************************/

#include <string>
#include <iostream>

#include <stdmix.h>
#include <mixio.h>
#include "qslim.h"

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

using namespace amrex;


#define SIZET INTEGER4
#define INTEGER4 int // from TECXXX.h



// Configuration variables and initial values
//
unsigned int face_target = 0;
bool will_use_fslim = false;
int placement_policy = MX_PLACE_OPTIMAL;
double boundary_weight = 1000.0;
int weighting_policy = MX_WEIGHT_AREA;
bool will_record_history = false;
double compactness_ratio = 0.0;
double meshing_penalty = 1.0;
bool will_join_only = false;
bool be_quiet = false;
OutputFormat output_format = MEF;
InputFormat input_format = iMEF;
char *output_filename = NULL;

// Globally visible variables
//
MxSMFReader *smf = NULL;
MxStdModel *m = NULL;
MxStdModel *m_orig = NULL;
MxQSlim *slim = NULL;
MxEdgeQSlim *eslim = NULL;
MxFaceQSlim *fslim = NULL;
QSlimLog *history = NULL;
MxDynBlock<MxEdge> *target_edges = NULL;

const char *slim_copyright_notice =
"Copyright (C) 1998-2002 Michael Garland.  See \"COPYING.txt\" for details.";

const char *slim_version_string = "2.1";

static ostream *output_stream = NULL;

static
bool qslim_smf_hook(char *op, int, char *argv[], MxStdModel& m)
{
    if( streq(op, "e") )
    {
	if( !target_edges )
	    target_edges = new MxDynBlock<MxEdge>(m.vert_count() * 3);

	MxEdge& e = target_edges->add();

	e.v1 = atoi(argv[0]) - 1;
	e.v2 = atoi(argv[1]) - 1;

	return true;
    }

    return false;
}

bool (*unparsed_hook)(char *, int, char*[], MxStdModel&) = qslim_smf_hook;

void slim_print_banner(ostream& out)
{
    out << "QSlim surface simplification software." << endl
	<< "Version " << slim_version_string << " "
	<< "[Built " << __DATE__ << "]." << endl
	<< slim_copyright_notice << endl;
}

void slim_init()
{
    if( !slim )
    {
	if( will_use_fslim )
	    slim = fslim = new MxFaceQSlim(*m);
	else
	    slim = eslim = new MxEdgeQSlim(*m);
    }
    else
    {
	if( will_use_fslim )
	    fslim = (MxFaceQSlim *)slim;
	else
	    eslim = (MxEdgeQSlim *)slim;
    }

    slim->placement_policy = placement_policy;
    slim->boundary_weight = boundary_weight;
    slim->weighting_policy = weighting_policy;
    slim->compactness_ratio = compactness_ratio;
    slim->meshing_penalty = meshing_penalty;
    slim->will_join_only = will_join_only;

    if( eslim && target_edges )
    {
	eslim->initialize(*target_edges, target_edges->length());
    }
    else
	slim->initialize();

    if( will_record_history )
    {
	if( !eslim )
	    mxmsg_signal(MXMSG_WARN,
			 "History only available for edge contractions.");
	else
	{
	    history = new QSlimLog(100);
	    eslim->contraction_callback = slim_history_callback;
	}
    }
}

#define CLEANUP(x)  if(x) { delete x; x=NULL; }

void slim_cleanup()
{
    CLEANUP(smf);
    CLEANUP(m);
    CLEANUP(slim);
    eslim = NULL;
    fslim = NULL;
    CLEANUP(history);
    CLEANUP(target_edges);
    if( output_stream != &cout )
    	CLEANUP(output_stream);
}

static
void setup_output()
{
    if( !output_stream )
    {
	if( output_filename )
            output_stream = new ofstream(output_filename);
        else
	    output_stream = &cout;
    }
}

bool select_output_format(const char *fmt)
{
    bool h = false;

    if     ( streq(fmt, "mmf") ) { output_format = MMF; h = true; }
    else if( streq(fmt, "pm")  ) { output_format = PM;  h = true; }
    else if( streq(fmt, "log") ) { output_format = LOG; h = true; }
    else if( streq(fmt, "smf") ) output_format = SMF;
    else if( streq(fmt, "iv")  ) output_format = IV;
    else if( streq(fmt, "vrml")) output_format = VRML;
    else if( streq(fmt, "mef")) output_format = MEF;
    else return false;

    if( h ) will_record_history = true;

    return true;
}

bool select_input_format(const char *fmt)
{
    if     ( streq(fmt, "mef") ) { input_format = iMEF;}
    else if( streq(fmt, "smf") ) { input_format = iSMF;}
    else return false;
    return true;
}

void output_preamble()
{
    if( output_format==MMF || output_format==LOG )
	output_current_model();
}

void output_current_model()
{
    setup_output();

    MxSMFWriter writer;
    writer.write(*output_stream, *m);
}

class FABdata
{
public:
    FABdata(SIZET i, SIZET n)
        {fab.resize(Box(IntVect::TheZeroVector(),
                        IntVect(AMREX_D_DECL(i-1,0,0))),n);}
    Real& operator[](SIZET i) {return fab.dataPtr()[i];}
    FArrayBox fab;
    SIZET boxSize;
};

void output_mef()
{
    std::string label = "decimated surface";
    std::string vars = "X Y Z";

    const SIZET nVert = m->vert_count();
    const int nComp = BL_SPACEDIM;

    FABdata vertData(nVert,nComp);
    SIZET cnt = 0;
    for(int i=0; i<nVert; i++)
    {
	if( m->vertex_is_valid(i) )
        {
            for (int j=0; j<nComp; ++j)
            {
                vertData[cnt+j] = m->vertex(i)[j];
            }
            cnt += nComp;
        }
    }

    // m->face_count apparently still includes a bunch of invalid faces, count valid ones
    SIZET nFaceTot = m->face_count();
    SIZET nFace = 0;
    for (int i=0; i<nFaceTot; ++i)
        if( m->face_is_valid(i) )
            nFace ++;

    const int vertsPerFace = 3;
    Vector<SIZET> faceData(nFace * vertsPerFace);
    cnt = 0;
    for(int i=0; i<nFaceTot; i++)
    {
	if( m->face_is_valid(i) )
        {
            for (int j=0; j<vertsPerFace; ++j)
            {
                faceData[cnt+j] = m->face(i)[j] + 1;
            }
            cnt += vertsPerFace;
        }
    }

    *output_stream << label << endl;
    *output_stream << vars << endl;
    *output_stream << nFace << " " << vertsPerFace << endl;
    vertData.fab.writeOn(*output_stream);
    output_stream->write((char*)faceData.dataPtr(),sizeof(SIZET)*faceData.size());
}

static
void cleanup_for_output()
{
    // First, mark stray vertices for removal
    //
    for(uint i=0; i<m->vert_count(); i++)
	if( m->vertex_is_valid(i) && m->neighbors(i).length() == 0 )
	    m->vertex_mark_invalid(i);

	// Compact vertex array so only valid vertices remain
    m->compact_vertices();
}

void output_final_model()
{
    setup_output();

    switch( output_format )
    {
    case MMF:
	output_regressive_mmf(*output_stream);
	break;

    case LOG:
	output_regressive_log(*output_stream);
	break;

    case PM:
	output_progressive_pm(*output_stream);
	break;

    case IV:
	cleanup_for_output();
	output_iv(*output_stream);
	break;

    case VRML:
	cleanup_for_output();
	output_vrml(*output_stream);
	break;

    case SMF:
	cleanup_for_output();
	output_current_model();
	break;

    case MEF:
	cleanup_for_output();
	output_mef();
	break;
    }


}

static std::string parseTitle(std::istream& is);
static std::vector<std::string> parseVarNames(std::istream& is);

void read_mef(std::istream* is)
{
    const std::string title = parseTitle(*is);
    const std::vector<std::string> names = parseVarNames(*is);
    const int nComp = names.size();

    SIZET nElts;
    int MYLEN;
    (*is) >> nElts;
    (*is) >> MYLEN;

    FArrayBox vertFab;
    vertFab.readFrom((*is));
    Real* nodeData = vertFab.dataPtr();
    int nPts = vertFab.box().numPts();

    for (int i=0; i<nPts; ++i)
    {
        int offset = i*nComp;
        m->add_vertex((float)nodeData[offset],
                      (float)nodeData[offset+1],
                      (float)nodeData[offset+2]);
    }
    vertFab.clear();

    Vector<SIZET> faceData(nElts*MYLEN,0);
    (*is).read((char*)faceData.dataPtr(),sizeof(SIZET)*faceData.size());

    for (int i=0; i<nElts; ++i)
    {
        int offset = i*MYLEN;
        m->add_face((unsigned int)faceData[offset] - 1,
                    (unsigned int)faceData[offset+1] - 1,
                    (unsigned int)faceData[offset+2] - 1);
    }
    faceData.clear();
}


void input_file(const char *filename)
{
    if( streq(filename, "-") )
    {
        if (input_format == iSMF)
        {
            smf->read(cin, m);
        }
        else
        {
            read_mef(&cin);
        }
    }
    else
    {
        ifstream in(filename);
        if( !in.good() )
            mxmsg_signal(MXMSG_FATAL, "Failed to open input file", filename);
        if (input_format == iSMF)
        {
            smf->read(in, m);
        }
        else
        {
            read_mef(&in);
        }
        in.close();
    }
}

static
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return Tokenize(line,std::string(", "));
}

static std::string parseTitle(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return line;
}

static
MxDynBlock<char*> files_to_include(2);

void defer_file_inclusion(char *filename)
{
    files_to_include.add(filename);
}

void include_deferred_files()
{
    for(uint i=0; i<files_to_include.length(); i++)
	input_file(files_to_include[i]);
}

void slim_history_callback(const MxPairContraction& conx, float cost)
{
    history->add(conx);
}
