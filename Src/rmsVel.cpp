#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_DataServices.H>
#include <AMReX_WritePlotFile.H>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    // Open first plotfile header and create an amrData object pointing into it
    int nPlotFiles = pp.countval("infiles");
    AMREX_ALWAYS_ASSERT(nPlotFiles>0);
    Vector<std::string> plotFileNames; pp.getarr("infiles",plotFileNames,0,nPlotFiles);

    int nVars(3);
    Vector<std::string> whichVar(nVars);
    whichVar[0] = "x_velocity";
    whichVar[1] = "y_velocity";
    whichVar[2] = "z_velocity";

    Vector<int>  destFills(nVars);
    for (int c=0; c<nVars; c++)
	destFills[c] = c;

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    Vector<DataServices *> dataServicesPtrVector(nPlotFiles);                                         // DataServices array for each plot
    Vector<AmrData *>      amrDataPtrVector(nPlotFiles);                                              // DataPtrVector for each plot
    Vector<Real>           time(nPlotFiles);
    Vector<Real>           urms(nPlotFiles);

    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
	if (ParallelDescriptor::IOProcessor())
	    std::cout << "Loading " << plotFileNames[iPlot] << std::endl;
	
	dataServicesPtrVector[iPlot] = new DataServices(plotFileNames[iPlot], fileType);               // Populate DataServices array
	
	if( ! dataServicesPtrVector[iPlot]->AmrDataOk())                                               // Check AmrData ok
	    DataServices::Dispatch(DataServices::ExitRequest, NULL);                                    // Exit if not
	
	amrDataPtrVector[iPlot] = &(dataServicesPtrVector[iPlot]->AmrDataRef());                        // Populate DataPtrVector
	
	time[iPlot] = amrDataPtrVector[iPlot]->Time();
	
    }

    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {
    
	int finestLevel = amrDataPtrVector[iPlot]->FinestLevel();    
	int inFinestLevel(-1);    pp.query("finestLevel",inFinestLevel);
	if (inFinestLevel>-1 && inFinestLevel<finestLevel) {
	    finestLevel = inFinestLevel;
            if (ParallelDescriptor::IOProcessor())
	        std::cout << "Finest level: " << finestLevel << std::endl;
	}

	Vector<Real> probLo=amrDataPtrVector[iPlot]->ProbLo();
	Vector<Real> probHi=amrDataPtrVector[iPlot]->ProbHi();
	const Real *dx = amrDataPtrVector[iPlot]->DxLevel()[finestLevel].dataPtr();
	Real dxyz = dx[0]*dx[1]*dx[2];

	int ngrow(0);
	MultiFab mf;
        const BoxArray& ba = amrDataPtrVector[iPlot]->boxArray(finestLevel);
        DistributionMapping dm(ba);
	mf.define(ba, dm, nVars, ngrow);

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "Processing " << iPlot << "/" << nPlotFiles << std::endl;
	amrDataPtrVector[iPlot]->FillVar(mf, finestLevel, whichVar, destFills);
	for (int n=0; n<nVars; n++)
	    amrDataPtrVector[iPlot]->FlushGrids(amrDataPtrVector[iPlot]->StateNumber(whichVar[n]));

	Real vol(0), uxb(0), uyb(0), uzb(0), ux2(0), uy2(0), uz2(0);
	for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
	    const FArrayBox &myFab = mf[ntmfi];
	    Vector<const Real *> varPtr(nVars);
	    for (int v=0; v<nVars; v++)
		varPtr[v] = myFab.dataPtr(v);
            const Box& vbx = ntmfi.validbox();
	    const int  *lo  = vbx.smallEnd().getVect();
	    const int  *hi  = vbx.bigEnd().getVect();

	    int ix = hi[0]-lo[0]+1;
	    int jx = hi[1]-lo[1]+1;
	    int kx = hi[2]-lo[2]+1;
	    for (int k=0; k<kx; k++) {
		Real z=probLo[2] + dx[2]*(0.5+(Real)(k+lo[2]));
		for (int j=0; j<jx; j++) {
		    Real y=probLo[1] + dx[1]*(0.5+(Real)(j+lo[1]));
		    for (int i=0; i<ix; i++) {
			Real x=probLo[0] + dx[0]*(0.5+(Real)(i+lo[0]));
			int cell = (k*jx+j)*ix+i;
			Real ux = varPtr[0][cell];
			Real uy = varPtr[1][cell];
			Real uz = varPtr[2][cell];
			vol += dxyz;
			uxb += ux*dxyz;
			uyb += uy*dxyz;
			uzb += uz*dxyz;
			ux2 += ux*ux*dxyz;
			uy2 += uy*uy*dxyz;
			uz2 += uz*uz*dxyz;
		    }
		}
	    }
	}
	ParallelDescriptor::ReduceRealSum(vol);
	ParallelDescriptor::ReduceRealSum(uxb);
	ParallelDescriptor::ReduceRealSum(uyb);
	ParallelDescriptor::ReduceRealSum(uzb);
	ParallelDescriptor::ReduceRealSum(ux2);
	ParallelDescriptor::ReduceRealSum(uy2);
	ParallelDescriptor::ReduceRealSum(uz2);
	uxb /= vol;	uyb /= vol;	uzb /= vol;
	ux2 /= vol;	uy2 /= vol;	uz2 /= vol;
	urms[iPlot] = sqrt( ( (ux2-uxb*uxb) + (uy2-uyb*uyb) + (uz2-uzb*uzb) ) / 3. );
    }
    if (ParallelDescriptor::IOProcessor())
        std::cout << "   ...done." << std::endl;

    if (ParallelDescriptor::IOProcessor()) {
	FILE *file = fopen("RmsVel.dat","w");
	for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
	    fprintf(file,"%e %e\n",time[iPlot],urms[iPlot]);
	fclose(file);
    }


  }
  Finalize();
  return 0;
}

