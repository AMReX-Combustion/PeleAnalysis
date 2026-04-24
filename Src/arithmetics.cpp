#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "This tool applies a simple arithmetic operation to two field "
       "variables.\n\n"
    << "Usage:\n"
    << "  " << argv[0]
    << " infile=FILE inVarAName=NAME inVarBName=NAME outVarName=NAME"
       " operator=OP [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE                     AMReX plotfile\n"
    << "  inVarAName=NAME                 operand A\n"
    << "  inVarBName=NAME                 operand B\n"
    << "  outVarName=NAME                 result variable name\n"
    << "  operator=OP                     operation: add, subtract, multiply, "
       "divide\n\n"

    << "Options:\n"
    << "  outfile=FILE                    output plotfile [DEF: infile_OP]\n"
    << "  coord=INT                       coordinate system (0=Cartesian) "
       "[DEF: 0]\n"
    << "  checkDivByZero=1|0              enable/disable divide-by-zero check "
       "[DEF: 1]\n"
    << "  -h, --help                      show this help message\n\n"

    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }
    ParmParse pp;

    if (pp.contains("verbose")) {
      AmrData::SetVerbose(true);
    }

    std::string infileName;
    pp.get("infile", infileName);

    std::string oper;
    pp.get("operator", oper);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      oper == "add" || oper == "subtract" || oper == "multiply" ||
        oper == "divide",
      "operator must be one of: add, subtract, multiply, divide");

    std::string outfileName(getFileRoot(infileName) + "_" + oper);
    pp.query("outfile", outfileName);

    int checkDivByZero = 1;
    if (oper == "divide") {
      pp.query("checkDivByZero", checkDivByZero);
    }

    std::string inVarAName;
    pp.get("inVarAName", inVarAName);
    std::string inVarBName;
    pp.get("inVarBName", inVarBName);
    std::string outVarName;
    pp.get("outVarName", outVarName);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(infileName, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(
        DataServices::ExitRequest, NULL); // does not return
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;

    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }
    Print() << "\n";

    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));

    int nCompIn = amrData.NComp();
    int nCompOut = nCompIn + 1;
    int outVar_id = nCompOut - 1;

    Vector<std::string> inNames = amrData.PlotVarNames();

    auto idA = std::find(inNames.begin(), inNames.end(), inVarAName);
    if (idA == inNames.end()) {
      Abort("Variable " + inVarAName + " not found in file " + infileName);
    }
    int inVarA_id = std::distance(inNames.begin(), idA);

    auto idB = std::find(inNames.begin(), inNames.end(), inVarBName);
    if (idB == inNames.end()) {
      Abort("Variable " + inVarBName + " not found in file " + infileName);
    }
    int inVarB_id = std::distance(inNames.begin(), idB);

    Vector<std::string> outNames = amrData.PlotVarNames();
    outNames.push_back(outVarName);

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    int coord = 0;
    pp.query("coord", coord);

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      outdata[lev].define(ba, dm, nCompOut, 0);

      Print() << "Reading data for level " << lev << std::endl;
      for (int i = 0; i < nCompIn; ++i) {
        outdata[lev].ParallelCopy(amrData.GetGrids(lev, i), 0, i, 1);
      }
      Print() << "Data has been read for level " << lev << std::endl;

      MultiFab::Copy(outdata[lev], outdata[lev], inVarA_id, outVar_id, 1, 0);
      if (oper == "add") {
        MultiFab::Add(outdata[lev], outdata[lev], inVarB_id, outVar_id, 1, 0);
      } else if (oper == "subtract") {
        MultiFab::Subtract(
          outdata[lev], outdata[lev], inVarB_id, outVar_id, 1, 0);
      } else if (oper == "multiply") {
        MultiFab::Multiply(
          outdata[lev], outdata[lev], inVarB_id, outVar_id, 1, 0);
      } else if (oper == "divide") {
        if (checkDivByZero) {
          ReduceOps<ReduceOpLogicalOr> reduce_op;
          ReduceData<int> reduce_data(reduce_op);
          using ReduceTuple = typename decltype(reduce_data)::Type;
          for (MFIter mfi(outdata[lev], TilingIfNotGPU()); mfi.isValid();
               ++mfi) {
            const Box& bx = mfi.tilebox();
            Array4<Real const> const arr = outdata[lev].const_array(mfi);
            reduce_op.eval(
              bx, reduce_data,
              [arr, inVarB_id] AMREX_GPU_HOST_DEVICE(int i, int j, int k)
                -> ReduceTuple { return {arr(i, j, k, inVarB_id) == 0.0}; });
          }
          int hasZero = amrex::get<0>(reduce_data.value());
          ParallelDescriptor::ReduceIntMax(hasZero);
          AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !hasZero, "Division by zero: " + inVarBName +
                        " contains zero values at level " +
                        std::to_string(lev));
        }
        MultiFab::Divide(
          outdata[lev], outdata[lev], inVarB_id, outVar_id, 1, 0);
      }
    }

    Print() << "Writing new data to " << outfileName << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      int r = amrData.RefRatio()[lev];
      refRatios[lev] = IntVect(AMREX_D_DECL(r, r, r));
    }
    amrex::WriteMultiLevelPlotfile(
      outfileName, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0,
      isteps, refRatios);
  }
  Finalize();
  return 0;
}
