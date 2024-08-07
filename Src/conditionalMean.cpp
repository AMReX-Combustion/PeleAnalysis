#include <string>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<filename>\n"
	      << "   binComp=i            : Variable id number to condition on\n"
	      << "   AvgComps=j k l       : Variable id numbers to average\n"
	      << "   min=m;  max=m        : min/max values for bins%i\n"
              << "   nBins=n              : Number of bins in PDF (default=64)\n"
	      << "   finestLevel=n        : Finest level at which to evaluate PDF\n"
	      << "   outSuffix=str        : Suffix to add to the pltfile name as an alt dir for results (default="")\n"
	      << "   aja=true/false       : Put the header in a separate file for gnuplot/matlab (default=false)\n"
              << "   infile= plt1 plt2    : List of plot files" << std::endl;
    exit(1);
}

typedef BaseFab<int> IntFab;
typedef FabArray<IntFab> IntMF;

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool isioproc = ParallelDescriptor::IOProcessor();
    int ioproc = ParallelDescriptor::IOProcessorNumber();
    int verbose=0; pp.query("verbose",verbose);
    if ( !(isioproc) )
        verbose=0;

    if (verbose>1)
        AmrData::SetVerbose(true);

    // Get plot file
    int nPlotFiles(pp.countval("infile"));
    if(nPlotFiles <= 0) {
        std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
        std::cerr << "Exiting." << std::endl;
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    if (verbose)
        std::cout << "Processing " << nPlotFiles << " plotfiles..." << std::endl;

    std::string outSuffix = ""; pp.query("outSuffix",outSuffix);

    // Make an array of srings containing paths of input plot files
    Vector<std::string> plotFileNames(nPlotFiles);
    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
        pp.get("infile", plotFileNames[iPlot], iPlot);
        if (verbose)
            std::cout << "   " << plotFileNames[iPlot] << std::endl;
    }

    // Get finest level argument
    int finestLevel(-1);
    pp.query("finestLevel",finestLevel);

    // Number of bins
    int nBins(64); pp.query("nBins",nBins);

    // Variables to bin, and to average
    int binComp = -1; pp.get("binComp",binComp);
    int nAvgComps = pp.countval("avgComps");
    Vector<int> avgComps;
    if (nAvgComps > 0)
    {
        avgComps.resize(nAvgComps);
        pp.getarr("avgComps",avgComps,0,nAvgComps);
    }
    else
    {
        amrex::Abort("need to specify avgComps");
    }
    Vector<Real>    binVals(nBins*nAvgComps,0);
    Vector<Real>  binValsSq(nBins*nAvgComps,0);
    bool writeBinMinMax = false; pp.query("writeBinMinMax",writeBinMinMax);
    Vector<Real> binMinVals;
    Vector<Real> binMaxVals;
    Vector<int>  binHits(nBins,0);
    if (writeBinMinMax)
    {
        binMinVals.resize(nBins*nAvgComps,0);
        binMaxVals.resize(nBins*nAvgComps,0);
    }
    
    Real binMin=0; pp.get("binMin",binMin);
    Real binMax=1; pp.get("binMax",binMax);
    if (binMax <= binMin)
        amrex::Abort("Bad bin min,max");


    bool floor=false; pp.query("floor",floor);
    bool ceiling=false; pp.query("ceiling",ceiling);
    
    Real domainVol = -1;
    Vector<int> weights;
    Vector<std::string> compNames;
    Vector<int> destFillComps(avgComps.size()+1);
    for (int i=0; i<destFillComps.size(); ++i)
        destFillComps[i] = i+1;// zero slot for mask
    int nComp = destFillComps.size() + 1;
    
    Vector<Real> dataIV(nComp);

    Vector<Real> bbll,bbur;
    if (int nx=pp.countval("bounds"))
    {
        Vector<Real> barr;
        pp.getarr("bounds",barr,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(barr.size()==2*d);
        bbll.resize(d);
        bbur.resize(d);
        for (int i=0; i<d; ++i)
        {
            bbll[i] = barr[i];
            bbur[i] = barr[d+i];
        }
    }

    bool aja(false);
    pp.query("aja",aja);
    if (aja && isioproc)
	std::cout << "Output for aja" << std::endl;

    // Loop over files
    for (int iPlot=0; iPlot<nPlotFiles; iPlot++)
    {
        // Open file and get an amrData pointer
        const std::string& infile = plotFileNames[iPlot];
        if (verbose)
            std::cout << "\nOpening " << infile << "..." << std::endl;
        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(infile, fileType);
        if (!dataServices.AmrDataOk())
            DataServices::Dispatch(DataServices::ExitRequest, NULL);

        AmrData& amrData = dataServices.AmrDataRef();
        int ngrow = 0;

        Box domain;
        if (iPlot==0)
        {
            compNames.resize(avgComps.size()+1);
            compNames[0] = amrData.PlotVarNames()[binComp];
            for (int i=0; i<avgComps.size(); ++i)
            {
                int comp = avgComps[i];
                if (comp<0 || comp>=amrData.NComp())
                    amrex::Abort("Bad comp: " + comp);
                compNames[i+1] = amrData.PlotVarNames()[comp];
            }

            domain = amrData.ProbDomain()[0];

            if (bbll.size()==BL_SPACEDIM  && bbur.size()==BL_SPACEDIM)
            {
                // Find coarse-grid coordinates of bounding box, round outwardly
                for (int i=0; i<BL_SPACEDIM; ++i) {
                    const Real dx = amrData.ProbSize()[i] / amrData.ProbDomain()[0].length(i);
                    domain.setSmall(i,
                                    std::max(domain.smallEnd()[i],
                                             (int)((bbll[i]-amrData.ProbLo()[i]+.0001*dx)/dx)));
                    domain.setBig(i,
                                  std::min(domain.bigEnd()[i],
                                           (int)((bbur[i]-amrData.ProbLo()[i]-.0001*dx)/dx)));
                }
            }
            domainVol = domain.d_numPts();

            if (finestLevel<0)
                finestLevel = amrData.FinestLevel();

            weights.resize(finestLevel+1,1);
            for (int i=finestLevel-1; i>=0; --i)
            {
                int rat = amrData.RefRatio()[i];
                weights[i] = weights[i+1];
                for (int d=0; d<BL_SPACEDIM; ++d)
                    weights[i] *= rat;
            }
        }
        
        // Build boxarrays for fillvar call
        Box levelDomain = domain;
        int thisFinestLevel = std::min(finestLevel,amrData.FinestLevel());
        int Nlevels = thisFinestLevel+1;
        Vector<BoxArray> bas(Nlevels);
        for (int iLevel=0; iLevel<=thisFinestLevel; ++iLevel)
        {
            BoxArray baThisLev = amrex::intersect(amrData.boxArray(iLevel),levelDomain);
            BoxList blThisLev(baThisLev);
            blThisLev.simplify();
            blThisLev.simplify();
            baThisLev = BoxArray(blThisLev);
            
            if (baThisLev.size() > 0) {
	        //bas.set(iLevel,baThisLev);
		bas[iLevel] = baThisLev;
		bas[iLevel].maxSize(32);
                if (iLevel < thisFinestLevel) {
                    levelDomain.refine(amrData.RefRatio()[iLevel]);
                }
            }
            else
            {
                bas.resize(iLevel);
            }
        }

        int ratio = 1;
        for (int iLevel=0; iLevel<=thisFinestLevel; ++iLevel)
        {
            if (bas.size()>iLevel)
            {
                DistributionMapping dmap(bas[iLevel]);
                MultiFab mf(bas[iLevel], dmap, nComp, ngrow);
                amrData.FillVar(mf, iLevel, compNames, destFillComps);
                mf.setVal(1,0,1); // initlaize mask to value

                // zero out covered data
                if (iLevel<thisFinestLevel)
                {
                    ratio = amrData.RefRatio()[iLevel];
                    BoxArray baf = BoxArray(bas[iLevel+1]).coarsen(ratio);
                    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
                    {
                        const Box& box = mfi.validbox();
                        FArrayBox& fab = mf[mfi];
                        std::vector< std::pair<int,Box> > isects = baf.intersections(box);
                        for (int ii = 0; ii < isects.size(); ii++)
                            fab.setVal(0,isects[ii].second,0,1);
                    }
                }

                for(MFIter mfi(mf); mfi.isValid(); ++mfi)
                {
                    FArrayBox& fab = mf[mfi];
                    const Box& box = mfi.validbox();
                    for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv))
                    {
                        fab.getVal(dataIV.dataPtr(),iv);
                        if (dataIV[0] != 0) // not covered
                        {
                            Real binVal = dataIV[1];
                            if (binVal>=binMin && binVal<binMax)
                            {
                                int myBin = (int)(nBins*(binVal-binMin)/(binMax-binMin));
                                if (myBin<0 || myBin>=binHits.size())
                                    amrex::Abort("Bad bin");
                                int myWeight = weights[iLevel];
                                for (int j=0; j<nAvgComps; ++j)
                                {
                                    Real val = dataIV[j+2];
                                    int binIdx = myBin*nAvgComps + j;
                                    binVals[binIdx]    += myWeight * val;
                                    binValsSq[binIdx]  += myWeight * val * val;
                                    if (writeBinMinMax)
                                    {
                                        if (binHits[myBin]==0) {
                                            binMinVals[binIdx] = val;
                                            binMaxVals[binIdx] = val;
                                        }
                                        binMinVals[binIdx] = std::min(val,binMinVals[binIdx]);
                                        binMaxVals[binIdx] = std::max(val,binMaxVals[binIdx]);
                                    }
                                }
                                binHits[myBin] += myWeight;
                            }
                        }
                    }
                } // MFI
            } // if
        } // Level
    } // iPlot
        
    ParallelDescriptor::ReduceIntSum(binHits.dataPtr(),binHits.size(),ioproc);
    ParallelDescriptor::ReduceRealSum(binVals.dataPtr(),binVals.size(),ioproc);
    ParallelDescriptor::ReduceRealSum(binValsSq.dataPtr(),binValsSq.size(),ioproc);
    if (writeBinMinMax)
    {
        ParallelDescriptor::ReduceRealMin(binMinVals.dataPtr(),binMinVals.size(),ioproc);
        ParallelDescriptor::ReduceRealMax(binMaxVals.dataPtr(),binMaxVals.size(),ioproc);
    }
    
    // Output result
    if (isioproc)
    {            
        std::string filename;
        // Output data for tecplot
        //filename = infile + outSuffix + "/CM_" + compNames[0] + ".dat";
	if (aja) {
	    filename = plotFileNames[0] + "/CM_" + compNames[0] + ".key";
	    std::cout << "Opening file " << filename << std::endl;
	} else {
	    // Default
	    filename = "CM_" + compNames[0] + ".dat";
	    std::cout << "Opening file " << filename << std::endl;
	}
        std::ofstream ofs(filename.c_str());
        std::string variables = "VARIABLES = " + compNames[0];
	for (int i=1; i<compNames.size(); ++i)
	    variables += " " + compNames[i] + "_sum";
	for (int i=1; i<compNames.size(); ++i)
	    variables += " " + compNames[i] + "_sumSq";
	for (int i=1; i<compNames.size(); ++i)
	    variables += " " + compNames[i] + "_avg";
	for (int i=1; i<compNames.size(); ++i)
	    variables += " " + compNames[i] + "_std";
        if (writeBinMinMax)
        {
            for (int i=1; i<compNames.size(); ++i)
                variables += " " + compNames[i] + "_min";
            for (int i=1; i<compNames.size(); ++i)
                variables += " " + compNames[i] + "_max";
        }
	variables += " N ";
	variables += " p ";
        variables += '\n';
        ofs << variables.c_str();
        ofs << "ZONE I=" << nBins << " DATAPACKING=POINT\n";
	if (aja) {
	    ofs.close();
	    filename = plotFileNames[0] + "/CM_" + compNames[0] + ".dat";
	    std::cout << "Opening file " << filename << std::endl;
	    ofs.open(filename.c_str());
	}
        const Real dv = (binMax - binMin)/nBins;
        int ntot = 0;
	for (int i=0; i<nBins; i++)
	{
	    ntot += binHits[i];
	}
        for (int i=0; i<nBins; i++) {
            const Real v = binMin + dv*(0.5+(Real)i);
            ofs << v << " ";
	    // Sum
            for (int j=0; j<nAvgComps; j++) {
                ofs << binVals[i*nAvgComps + j] << " ";
            }
	    // SumSq
            for (int j=0; j<nAvgComps; j++) {
                ofs << binValsSq[i*nAvgComps + j] << " ";
            }
	    if (binHits[i]>0) {
		// Avg
		for (int j=0; j<nAvgComps; j++) {
		    ofs << binVals[i*nAvgComps + j]/(Real)binHits[i] << " ";
		}
		// Std Dev
		for (int j=0; j<nAvgComps; j++) {
		    int idx = i*nAvgComps + j;
		    Real bh = (Real)binHits[i];
                    ofs << std::sqrt( (binValsSq[idx]/bh) - (binVals[idx]/bh)*(binVals[idx]/bh) ) << " ";
		}
	    } else {
		// Avg & Std Dev
		for (int j=0; j<nAvgComps*2; j++) {
                    ofs << "0.0 ";
		}
	    }
	    if (writeBinMinMax)
            {
                for (int j=0; j<nAvgComps; j++) {
                    ofs << binMinVals[i*nAvgComps + j] << " ";
                }
                for (int j=0; j<nAvgComps; j++) {
                    ofs << binMaxVals[i*nAvgComps + j] << " ";
                }
            }
	    ofs << Real(binHits[i]) << " ";
            ofs << Real(binHits[i])/ntot << '\n';
        }
        cout << "total bins: " << ntot << endl;
        ofs.close();
        
    } // IOProcessor
    
    amrex::Finalize();
    return 0;
}

