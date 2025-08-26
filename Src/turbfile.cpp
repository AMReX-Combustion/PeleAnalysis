#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

#if AMREX_SPACEDIM == 2
extern "C"
{
    void FLIPROWSY(const Real* dat, ARLIM_P(lo),ARLIM_P(hi));
};
#endif

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " ifile=<pltfile> ofile=<turbname>\n";
    exit(1);
}

static
void
Extend (FArrayBox& xfab,
        FArrayBox& vfab,
        const Box& dm)
{
    Box tbx = vfab.box();

    tbx.setBig(0, dm.bigEnd(0) + 3);

    const int ygrow = AMREX_SPACEDIM==3 ? 3 : 1;

    tbx.setBig(1, dm.bigEnd(1) + ygrow);

    xfab.resize(tbx,1);

    xfab.copy(vfab);
    vfab.shift(0, dm.length(0));
    xfab.copy(vfab);
    vfab.shift(1, dm.length(1));
    xfab.copy(vfab);
    vfab.shift(0, -dm.length(0));
    xfab.copy(vfab);
}

int
main (int   argc,
      char* argv[])
{
    Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    if(ParallelDescriptor::NProcs() > 1) {
      Abort("Not suitable for parallel");
    }
    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    std::string ifile;
    pp.get("ifile",ifile);

    std::string ofile;
    pp.get("ofile",ofile);

    std::cout << "Reading " << ifile << " ... " << std::flush;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(ifile, fileType);
    std::cout << "done" << std::endl;

    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = dataServices.AmrDataRef();

    std::string TurbDir = ofile;

    if (ParallelDescriptor::IOProcessor())
        if (!UtilCreateDirectory(TurbDir, 0755))
            CreateDirectoryFailed(TurbDir);

    std::string Hdr = TurbDir; Hdr += "/HDR";
    std::string Dat = TurbDir; Dat += "/DAT";

    std::ofstream ifsd, ifsh;

    ifsh.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsh.good())
        FileOpenFailed(Hdr);

    ifsd.open(Dat.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsd.good())
        FileOpenFailed(Dat);

    int dir         = AMREX_SPACEDIM - 1; pp.query("dir",dir);
    int dir_arr[3];
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      d == dir ? dir_arr[d] = 0 : dir_arr[d] = 1;      
    }
    
    const int          finestLevel = amrData.FinestLevel();
    const Box          dm          = amrData.ProbDomain()[finestLevel];
    Box                xdm         = dm;
    const Vector<Real> dx          = amrData.DxLevel()[finestLevel];
    IntVect            sm          = dm.smallEnd();
    IntVect            bg          = dm.bigEnd();
    std::string        names[3]    = { "x_velocity", "y_velocity", "z_velocity" };
    FArrayBox          xfab;
    //
    // Write the first part of the header.
    // Note that this is solely for periodic style inflow files.
    //
#if AMREX_SPACEDIM==2
    xdm.setBig(0, dm.bigEnd(0) + 3);
    xdm.setBig(1, dm.bigEnd(1) + 1);

    ifsh << xdm.length(0) << ' '
         << xdm.length(1) << ' '
         << 1             << '\n';

    ifsh << amrData.ProbSize()[0] + 2*dx[0] << ' '
         << amrData.ProbSize()[1]           << ' '
         << 0                               << '\n';

    ifsh << 1 << ' ' << 1 << ' ' << 0 << '\n';
#elif AMREX_SPACEDIM==3
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      xdm.setBig(d, dm.bigEnd(d) + 1 + 2*dir_arr[d]);
    }

    ifsh << xdm.length(0) << ' '
         << xdm.length(1) << ' '
         << xdm.length(2) << '\n';

    ifsh << amrData.ProbSize()[0] + 2*dir_arr[0]*dx[0] << ' '
         << amrData.ProbSize()[1] + 2*dir_arr[1]*dx[1] << ' '
         << amrData.ProbSize()[2] + 2*dir_arr[2]*dx[2] << '\n';

    ifsh << 1 << ' ' << 1 << ' ' << 1 << '\n';
#else
    Abort("Only 2-D & 3-D supported");
#endif

    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        std::cout << "Loading component " << d << " ... " << std::flush;

        BoxArray ba(1);

#if AMREX_SPACEDIM==3
        //
        // In 3-D we work on one cell wide Z-planes.
        // We first do the lo AMREX_SPACEDIM plane.
        // And then all the other planes in xhi -> xlo order.
        //
        // In 2-D we only write a single x-y plane of data per component.
        //
        bg[dir] = sm[dir];
#endif
        Box bx(sm,bg);
        ba.set(0,bx);
        MultiFab TMP(ba,DistributionMapping(ba),1,0);

        amrData.FillVar(TMP, amrData.FinestLevel(), names[d], 0);
        amrData.FlushGrids(amrData.StateNumber(names[d]));

        Extend(xfab, TMP[0], dm);
        //
        // Write current position of data file to header file.
        //
        ifsh << ifsd.tellp() << std::endl;
#if AMREX_SPACEDIM==2
        //
        // Write the FAB to the data file after flipping rows in Y direction.
        //
        FLIPROWSY(xfab.dataPtr(), ARLIM(xfab.loVect()), ARLIM(xfab.hiVect()));
        xfab.writeOn(ifsd);
#elif AMREX_SPACEDIM==3

        xfab.writeOn(ifsd);

        std::cout << "xfab min,max: " << xfab.min() << ' ' << xfab.max() << ' ' << bx << '\n';
        //
        // Now do all the planes in dirhi -> dirlo order.
        //
        for (int i = dm.bigEnd(dir); i >= dm.smallEnd(dir); i--)
        {
            sm[dir] = i;
            bg[dir] = i;
            Box bx(sm,bg);
            ba.set(0,bx);
            MultiFab TMP(ba,DistributionMapping(ba),1,0);

            amrData.FillVar(TMP, amrData.FinestLevel(), names[d], 0);
            amrData.FlushGrids(amrData.StateNumber(names[d]));

            Extend(xfab, TMP[0], dm);
	    
            //
            // Write current position of data file to header file.
            //
            ifsh << ifsd.tellp() << std::endl;
            //
            // Write the FAB to the data file.
            //
            xfab.writeOn(ifsd);

            std::cout << "xfab i,min,max: " << i<< ' ' << xfab.min() << ' ' << xfab.max() << ' ' << bx << '\n';
        }
#endif
        std::cout << "done" << std::endl;
    }

    Finalize();
}
