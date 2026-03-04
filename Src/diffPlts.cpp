// takes two plotfiles and calculates the difference

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " infile1=FILE infile2=FILE [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile1=FILE       First AMReX plotfile\n"
            << "  infile2=FILE       Second AMReX plotfile\n\n"

            << "Options:\n"
            << "  -h, --help         Show this help message\n\n"
            << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
            << "the documentation.\n";

  std::exit(1);
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  if (argc < 2) {
    print_usage(argc, argv);
  } else if (
    (std::strcmp(argv[1], "-h") == 0) ||
    (std::strcmp(argv[1], "--help") == 0)) {
    print_usage(argc, argv);
  }

  // get infile names and count
  ParmParse pp;
  int nfiles(pp.countval("infiles"));
  Vector<std::string> infiles(nfiles);
  pp.getarr("infiles", infiles);

  if (nfiles != 2) {
    amrex::Abort("Tool is only designed for two infiles, please adjust");
  }
  // get and count variables to copy
  int nvars(pp.countval("vars"));
  Vector<std::string> vars(nvars);
  pp.getarr("vars", vars);
  Vector<std::string> new_vars = vars;
  Vector<std::string> names;
  // outfile name
  std::string outfile = infiles[0] + "_diff";
  pp.query("outfile", outfile);

  std::string diff_type = "absolute";
  pp.query("diff_type", diff_type);

  AMREX_ALWAYS_ASSERT((diff_type == "absolute") || (diff_type == "relative"));

  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

  // setting up for reading pltfiles
  DataServices dataServices0(infiles[0], fileType);
  DataServices dataServices1(infiles[1], fileType);

  if (!dataServices0.AmrDataOk() || !dataServices1.AmrDataOk())
    DataServices::Dispatch(DataServices::ExitRequest, NULL);

  AmrData& amrData0 = dataServices0.AmrDataRef();
  AmrData& amrData1 = dataServices1.AmrDataRef();

  // getting the finest level (or whatever user sets)
  // aborting if finest level dont match
  int finestLevel = std::min(amrData0.FinestLevel(), amrData1.FinestLevel());
  finestLevel = amrData0.FinestLevel();
  pp.query("finestLevel", finestLevel);

  if (
    (finestLevel > amrData1.FinestLevel()) ||
    (finestLevel > amrData0.FinestLevel()))
    amrex::Abort(
      "Requested finest level exeeds finest level of one of the plotfiles");

  int Nlev = finestLevel + 1;

  // setting up periodicity
  Vector<int> is_per(AMREX_SPACEDIM, 0);
  pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
  Print() << "Periodicity assumed for this case: ";
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    Print() << is_per[idim] << " ";
  }
  Print() << "\n";

  // setting up pltfile boxes
  RealBox rb(&(amrData0.ProbLo()[0]), &(amrData0.ProbHi()[0]));
  Vector<MultiFab*> fileData0(Nlev);
  Vector<MultiFab*> fileData1(Nlev);
  Vector<Geometry> geoms(Nlev);
  int coord = 0;
  // create the geometry of MultiFabs, the "boxArray" Parent does'nt matter
  for (int lev = 0; lev < Nlev; ++lev) {
    const DistributionMapping dm(amrData0.boxArray(lev));
    geoms[lev] = Geometry(amrData0.ProbDomain()[lev], &rb, coord, &(is_per[0]));
    fileData0[lev] = new MultiFab(amrData0.boxArray(lev), dm, nvars, 0);
    fileData1[lev] = new MultiFab(amrData1.boxArray(lev), dm, nvars, 0);
  }

  // data structure to create identical var tables
  const Vector<std::string>& plot0_VarNames = amrData0.PlotVarNames();
  const Vector<std::string>& plot1_VarNames = amrData1.PlotVarNames();

  Vector<int> idcomp0;
  Vector<int> idcomp1;

  // loop through all available vars
  for (int i = 0; i < nvars; i++) {
    const std::string& var = vars[i];
    int idx0 = -1;
    int idx1 = -1;

    // loop through first infile to find var number i
    for (int j = 0; j < plot0_VarNames.size(); j++) {
      // if var was found, save the position in idx
      if (plot0_VarNames[j] == var) {
        idx0 = j;
        break;
      }
    }

    // loop through second infile to finde var number i
    for (int k = 0; k < plot1_VarNames.size(); k++) {
      // if var was found, save the position in idx
      if (plot1_VarNames[k] == var) {
        idx1 = k;
        break;
      }
    }

    // if first idx is unchanged: var was not found, abort
    if (idx0 < 0) {
      Print() << "Error: Variable " << var << " not in first infile \n";
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // if second idx is unchanged: var was not found, abort
    if (idx1 < 0) {
      Print() << "Error: Variable " << var << " not in second infile \n";
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // save the position for each file, in which the var lays
    idcomp0.push_back(idx0);
    idcomp1.push_back(idx1);
    // now we have two tables with indexes ponting to the right var for each
    // infile
  }

  // fill in the Multifabs with the infile Data
  for (int lev = 0; lev < Nlev; lev++) {
    for (int i = 0; i < idcomp0.size(); i++) {
      fileData0[lev]->ParallelCopy(amrData0.GetGrids(lev, idcomp0[i]), 0, i, 1);
      fileData1[lev]->ParallelCopy(amrData1.GetGrids(lev, idcomp1[i]), 0, i, 1);
    }
  }

  // creating an array containing the var names

  std::string suffix;

  if (diff_type == "absolute") {
    suffix = "_diff";
  } else if (diff_type == "relative") {
    suffix = "_rel_diff";
  }

  for (int i = 0; i < idcomp0.size(); i++) {
    names.push_back(plot0_VarNames[idcomp0[i]] + suffix);
  }

  // calculate the difference
  for (int lev = 0; lev < Nlev; ++lev) {
    fileData1[lev]->Subtract(*fileData1[lev], *fileData0[lev], 0, 0, nvars, 0);
    if (diff_type == "relative") {
      fileData1[lev]->Divide(*fileData1[lev], *fileData0[lev], 0, 0, nvars, 0);
    }
  }

  // write pltfile
  Vector<int> isteps(Nlev, 0);
  Vector<IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
  amrex::WriteMultiLevelPlotfile(
    outfile, Nlev, GetVecOfConstPtrs(fileData1), names, geoms, 0.0, isteps,
    refRatios);

  amrex::Finalize();
  return 0;
}
