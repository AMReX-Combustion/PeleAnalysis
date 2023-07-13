
#include <AMReX_ParmParse.H>
#include <AMReX_AmrData.H>
#include <AMReX_DataServices.H>

using namespace amrex;

Real VolWgtAvgC(const FArrayBox& Val, int valComp,
                const FArrayBox& Vol, int volComp)
{
  Box box = Val.box() & Vol.box();
  FArrayBox tmp(box,1);
  tmp.copy<RunOn::Device>(Val,box,valComp,box,0,1);
  tmp.mult<RunOn::Device>(Vol,box,volComp,0,1);
  return tmp.sum<RunOn::Device>(box,0,1) / Vol.sum<RunOn::Device>(box,volComp,1);
}

Real VolWgtAvg(const FArrayBox& Val,
               const FArrayBox& Vol)
{
  Box box = Val.box() & Vol.box();
  FArrayBox tmp(box,1);
  tmp.copy<RunOn::Device>(Val);
  tmp.mult<RunOn::Device>(Vol,0,0);
  return tmp.sum<RunOn::Device>(0) / Vol.sum<RunOn::Device>(0);
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    ParmParse pp;

    if (pp.contains("help")) {
      //print_usage(argc,argv);
    }

    std::string infile;
    pp.get("infile",infile);

    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    AMREX_ALWAYS_ASSERT(finestLevel >= 0 && finestLevel<=amrData.FinestLevel());

    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));

    int nc = pp.countval("comps");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nc==6,"comps must be a list of 6 integers (adv_0, adv_1, Var1, Var2, Var3, vfrac");
    Vector<int> comps(6);
    pp.getarr("comps",comps,0,nc);

    Box subbox = amrData.ProbDomain()[finestLevel];
    Vector<int> inBox;

    if (int nx=pp.countval("box"))
    {
        pp.getarr("box",inBox,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(inBox.size()==2*d);
        subbox=Box(IntVect(D_DECL(inBox[0],inBox[1],inBox[2])),
                   IntVect(D_DECL(inBox[d],inBox[d+1],inBox[d+2])),
                   IndexType::TheCellType());
    }

    subbox &= amrData.ProbDomain()[finestLevel];
    
    Vector<std::string> names(comps.size());
    Vector<int> fillComps(comps.size());
    
    for (int i=0; i<comps.size(); ++i)
    {
      names[i] = amrData.PlotVarNames()[comps[i]];
      fillComps[i] = i;
    }

    int planeCoord = 0; pp.get("planeCoord",planeCoord);
    const Vector<Real>& plo = amrData.ProbLo();
    Vector<Real> dx = amrData.CellSize(finestLevel);
        
    IntVect se = subbox.smallEnd();
    IntVect be = subbox.bigEnd();
    int clo = se[planeCoord];
    int chi = be[planeCoord];

    Real v0min = 0;
    Real v0max = 1;
    int nBin0 = 64; pp.query("nBins",nBin0);
    int nBin1 = nBin0;
    Real v1min = 0;
    Real v1max = 1;
    Box binbox(IntVect(AMREX_D_DECL(0,0,0)),
               IntVect(AMREX_D_DECL(nBin0-1,nBin1-1,0)));
    int nBinR = nBin0;
    Box rbinbox(IntVect(AMREX_D_DECL(0,0,0)),
                IntVect(AMREX_D_DECL(0,nBinR-1,0)));

    FArrayBox outbins(binbox,2);
    FArrayBox outrbins(rbinbox,2);
    outbins.setVal<RunOn::Device>(0);
    outrbins.setVal<RunOn::Device>(0);
    int nBinPlanes = 10; pp.query("nBinPlanes",nBinPlanes);
    int numPlanesSoFar = 0;
    
    GpuArray<Real,AMREX_SPACEDIM-1> na;
    int icnt = 0;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
      if (i!=planeCoord) {
        na[icnt++] = i;
      }
    }
    GpuArray<Real,AMREX_SPACEDIM> ploa = {AMREX_D_DECL(plo[0],plo[1],plo[2])};
    GpuArray<Real,AMREX_SPACEDIM> dxa  = {AMREX_D_DECL( dx[0], dx[1], dx[2])};
    Real R = amrData.ProbHi()[na[0]];

    // Only the I/O processor makes the directory if it doesn't already exist.
    std::string output_dir = "Output";
    if (ParallelDescriptor::IOProcessor())
      if (!amrex::UtilCreateDirectory(output_dir, 0755))
        amrex::CreateDirectoryFailed(output_dir);

    // Force other processors to wait till directory is built.
    ParallelDescriptor::Barrier();


    std::string out_mean_file(output_dir + "/" + "mean.dat"); pp.query("out_mean_file",out_mean_file);
    std::string out_mcmt_file(output_dir + "/" + "mcmt.dat"); pp.query("out_mcmt_file",out_mcmt_file);
    std::ofstream os_mean(out_mean_file.c_str());
    std::ofstream os_mcmt(out_mcmt_file.c_str());

    int nPlanesPerPass = nBinPlanes;
    bool parallelBin = 1; pp.query("parallelBin",parallelBin);

    if (parallelBin)
    {
      BoxList bl;
      for (int c=clo; c<=chi; c+=nPlanesPerPass)
      {
        int cbhi = amrex::min(c + nPlanesPerPass - 1, chi);
        se[planeCoord] = c;
        be[planeCoord] = cbhi;
        bl.push_back(Box(se,be));
      }
      BoxArray ba(bl);
      MultiFab data(ba,DistributionMapping(ba),comps.size(),0);
      int nBoxes = ba.size();

      Print() << "Reading data..." << std::endl;
      amrData.FillVar(data,finestLevel,names,fillComps);
      Print() << "Data has been read.  Processing..." << std::endl;

      Vector<Real> volWgtAvg2(nBoxes,0);
      Vector<Real> volWgtAvg3(nBoxes,0);
      Vector<Real> volWgtAvg4(nBoxes,0);
      Vector<Real> mcmt(nBoxes,0);
      Vector<Real> planeLoc(nBoxes,0);

      for (MFIter mfi(data); mfi.isValid(); ++mfi)
      {
        const Box& b = mfi.validbox();
        int idx = mfi.index();
        FArrayBox bins(binbox,2);   // Bins for adv0,adv1
        FArrayBox rbins(rbinbox,2); // Bins for radius
        Elixir e_bins = bins.elixir();   // Keep the fabs around until gpus done
        Elixir e_rbins = rbins.elixir(); // --ditto

        bins.setVal<RunOn::Device>(0);
        rbins.setVal<RunOn::Device>(0);

        const auto& f0a = data[mfi].array(0); // Array4 for adv0
        const auto& f1a = data[mfi].array(1); // Array4 for adv1
        const auto& vala = data[mfi].array(2); // Qty to bin in (adv0,adv1)
        const auto& vola = data[mfi].array(comps.size()-1); //Volume (=1 for full cells)
        const auto& outa = bins.array();   // Array4 for adv0,adv1 bins
        const auto& outra = rbins.array(); // Array4 for radial bins

        ParallelFor(b,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                      {
                        int bin0 = amrex::max(0,amrex::min(nBin0-1,(int)((nBin0-1)*(f0a(i,j,k) - v0min)/(v0max - v0min))));
                        int bin1 = amrex::max(0,amrex::min(nBin1-1,(int)((nBin0-1)*(f1a(i,j,k) - v1min)/(v1max - v1min))));
                        HostDevice::Atomic::Add(&outa(bin0,bin1,0,0),vala(i,j,k) * vola(i,j,k));
                        HostDevice::Atomic::Add(&outa(bin0,bin1,0,1),vola(i,j,k));

                        Real y = ploa[na[0]] + (j+0.5)*dxa[na[0]];
                        Real z = ploa[na[1]] + (k+0.5)*dxa[na[1]];
                        Real r = sqrt(y*y + z*z);
                        int binr = amrex::max(0,amrex::min(nBinR,(int)((nBinR-1)*(r/R))));
                        HostDevice::Atomic::Add(&outra(0,binr,0),vala(i,j,k)*vola(i,j,k));
                        HostDevice::Atomic::Add(&outra(0,binr,1),vola(i,j,k));
                      });

        amrex::Gpu::synchronize();

        ParallelFor(binbox,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                      {
                        if (outa(i,j,k,1) > 0) {
                          outa(i,j,k,0) = outa(i,j,k,0) / outa(i,j,k,1);
                        }
                      });

        ParallelFor(rbinbox,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                      {
                        if (outra(i,j,k,1) > 0) {
                          outra(i,j,k,0) = outra(i,j,k,0) / outra(i,j,k,1);
                        }
                      });

        amrex::Gpu::synchronize();
      
        planeLoc[idx] = plo[planeCoord] + ( 0.5*(b.smallEnd(planeCoord)+b.bigEnd(planeCoord))+0.5)*dx[planeCoord];

        volWgtAvg2[idx] = VolWgtAvgC(data[mfi], 2, data[mfi], comps.size()-1);
        volWgtAvg3[idx] = VolWgtAvgC(data[mfi], 3, data[mfi], comps.size()-1);
        volWgtAvg4[idx] = VolWgtAvgC(data[mfi], 4, data[mfi], comps.size()-1);
        mcmt[idx] = rbins.max<RunOn::Device>(0);
      }
      ParallelDescriptor::ReduceRealSum(volWgtAvg2.dataPtr(), nBoxes);
      ParallelDescriptor::ReduceRealSum(volWgtAvg3.dataPtr(), nBoxes);
      ParallelDescriptor::ReduceRealSum(volWgtAvg4.dataPtr(), nBoxes);
      ParallelDescriptor::ReduceRealSum(mcmt.dataPtr(),       nBoxes);
      ParallelDescriptor::ReduceRealSum(planeLoc.dataPtr(),   nBoxes);

      if (ParallelDescriptor::IOProcessor())
      {
        for (int i=0; i<nBoxes; ++i) {
          os_mcmt << planeLoc[i] << " " << mcmt[i] << std::endl;
          os_mean << planeLoc[i] << " " << volWgtAvg2[i] <<  " " << volWgtAvg3[i] <<  " " << volWgtAvg4[i] << std::endl;
        }
      }
    }
    else
    {
      for (int c=clo; c<=chi; c+=nPlanesPerPass) {

        int cbhi = amrex::min(c + nPlanesPerPass - 1, chi);
        Real planeLoc = plo[planeCoord] + ( 0.5*(c+cbhi)+0.5)*dx[planeCoord];
        Print() << "Planes " << c << ":" << cbhi << " at " << planeLoc << std::endl;
      
        se[planeCoord] = c;
        be[planeCoord] = cbhi;
        Box fabbox(se,be);
        Array<std::unique_ptr<FArrayBox>, 6> datav;
        if (ParallelDescriptor::IOProcessor()) {
          for (int ic=0; ic<comps.size(); ic++) {
            datav[ic] = std::unique_ptr<FArrayBox>(new FArrayBox(fabbox,1));
          }
        }

        // get the data
        for (int ic=0; ic<comps.size(); ++ic) {
          amrData.FillVar(datav[ic].get(),fabbox,finestLevel,names[ic],ParallelDescriptor::IOProcessorNumber());
        }

        if (ParallelDescriptor::IOProcessor()) {

          // for each cell, find the 2D bin, and increment the volume and var value in the result fab
          const auto& f0a = datav[0]->array();
          const auto& f1a = datav[1]->array();
          const auto& vala = datav[2]->array();
          const auto& vola = datav[comps.size()-1]->array();
          const auto& outa = outbins.array();
          const auto& outra = outrbins.array();

          ParallelFor(fabbox,
                      [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                        {
                          int bin0 = amrex::max(0,amrex::min(nBin0-1,(int)((nBin0-1)*(f0a(i,j,k) - v0min)/(v0max - v0min))));
                          int bin1 = amrex::max(0,amrex::min(nBin1-1,(int)((nBin0-1)*(f1a(i,j,k) - v1min)/(v1max - v1min))));
                          HostDevice::Atomic::Add(&outa(bin0,bin1,0,0),vala(i,j,k) * vola(i,j,k));
                          HostDevice::Atomic::Add(&outa(bin0,bin1,0,1),vola(i,j,k));

                          Real y = ploa[na[0]] + (j+0.5)*dxa[na[0]];
                          Real z = ploa[na[1]] + (k+0.5)*dxa[na[1]];
                          Real r = sqrt(y*y + z*z);
                          int binr = amrex::max(0,amrex::min(nBinR,(int)((nBinR-1)*(r/R))));
                          HostDevice::Atomic::Add(&outra(0,binr,0),vala(i,j,k)*vola(i,j,k));
                          HostDevice::Atomic::Add(&outra(0,binr,1),vola(i,j,k));
                        });
          amrex::Gpu::synchronize();
        
          numPlanesSoFar+=nPlanesPerPass;
          if (numPlanesSoFar>=nBinPlanes) {
          
            ParallelFor(binbox,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                          {
                            if (outa(i,j,k,1) > 0) {
                              outa(i,j,k,0) = outa(i,j,k,0) / outa(i,j,k,1);
                            }
                          });

            ParallelFor(rbinbox,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
                          {
                            if (outra(i,j,k,1) > 0) {
                              outra(i,j,k,0) = outra(i,j,k,0) / outra(i,j,k,1);
                            }
                          });

            amrex::Gpu::synchronize();

            std::string plane_file = Concatenate("Output/bins_",c,5) + ".fab";
            std::ofstream osf;
            osf.open(plane_file.c_str());
            outbins.writeOn(osf);
            osf.close();

            std::string rplane_file = Concatenate("Output/rbins_",c,5) + ".fab";
            std::ofstream osfr;
            osfr.open(rplane_file.c_str());
            outrbins.writeOn(osfr);
            osfr.close();

            os_mcmt << planeLoc << " " << outrbins.max<RunOn::Device>(0) << std::endl;

            numPlanesSoFar = 0;
            outbins.setVal<RunOn::Device>(0);
            outrbins.setVal<RunOn::Device>(0);
          }

          Real volWgtAvg2 = VolWgtAvg(*datav[2],*datav[comps.size()-1]);
          Real volWgtAvg3 = VolWgtAvg(*datav[3],*datav[comps.size()-1]);
          Real volWgtAvg4 = VolWgtAvg(*datav[4],*datav[comps.size()-1]);

          os_mean << planeLoc << " " << volWgtAvg2 <<  " " << volWgtAvg3 <<  " " << volWgtAvg4 << std::endl;

        }
      }
    }
    
    os_mean.close();
    os_mcmt.close();
    amrex::Finalize();
}
