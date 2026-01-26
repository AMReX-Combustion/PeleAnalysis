#include <string>
#include <iostream>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

#include <streamBinTubeStats.H>


static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " inputs infile=<s> isoCompName=<s> isoVal=<v> [options] \n\tOptions:\n";
  std::cerr << "\t     infile=<s> where <s> is a stream.bin file\n";
  std::cerr << "\t     writeStreamsToMatlab\n";
  std::cerr << "\t     writeTecplotSurfaceFromStream\n";
  std::cerr << "\t     writeBasic\n";
  std::cerr << "\t     avgComps: list of variables to average\n";
  std::cerr << "\t     intComps: list of variables to integrate\n";
  std::cerr << "\t     derComps: derived values. Can be:\n";
  std::cerr << "\t     --> flameThickness\n";
  std::cerr << "\t            reacTemp\n";
  std::cerr << "\t            prodTemp\n";
  std::cerr << "\t            tempGradVar\n";
  std::cerr << "\t     --> principalCurvatureZone\n";
  std::cerr << "\t            pkzLength\n";
  std::cerr << "\t            pkzMkVar\n";
  std::cerr << "\t            pkzGkVar\n";
exit(1);
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  if (ParallelDescriptor::NProcs()>1)
    Abort("Code is not yet parallel safe");
  
  if (argc < 2) {
    print_usage(argc,argv);
  }

  ParmParse pp;

  if (pp.contains("help")) {
    print_usage(argc,argv);
  }
  
  // read infile name from inputs
  std::string infile;
  pp.get("infile",infile);

  // parse what files to write
  int writeStreamsToMatlab(0);
  pp.query("writeStreamsToMatlab",writeStreamsToMatlab);
  int writeTecplotSurfaceFromStream(0);
  pp.query("writeTecplotSurfaceFromStream",writeTecplotSurfaceFromStream);
  int writeBasic(0);
  pp.query("writeBasic",writeBasic);
  int writeSurface= 1;
  pp.query("writeSurface", writeSurface);
  
  std::string fuelName="H2";
  pp.query("fuelName",fuelName);
  // declare size and data holders
  int nStreams, nElts, nPtsOnStream, nComps;
  Vector<std::string> variableNames;
  Vector<int> faceData;
  Vector<Vector<Real>> streamData;

  // Read file
  readStreamBin(infile, nStreams, nElts, nPtsOnStream, nComps,
		variableNames, faceData, streamData);
  // report
  Print() << "nStreams      = " << nStreams << std::endl;
  Print() << "nElements     = " << nElts << std::endl;
  Print() << "nPtsOnStreams = " << nPtsOnStream << std::endl;
  Print() << "nComps        = " << nComps << std::endl;
  Print() << "Variable names:" << std::endl;

  for (int iComp=0; iComp<nComps; iComp++) 
    Print() << "   " << iComp << ": " << variableNames[iComp] << std::endl;
  
  // let's write a matlab file for each variable
  if (writeStreamsToMatlab) {
    Print() << "Writing streams as matlab files..." << std::endl;
    writeStreamsMatlab(infile,nStreams,nPtsOnStream,nComps,variableNames,streamData);
  }
  
  // let's output a surface
  if (writeTecplotSurfaceFromStream) {
    Print() << "Writing a tecplot surface..." << std::endl;
    writeSurfaceFromStreamTecplot(infile,nStreams,nElts,nPtsOnStream,nComps,
				  variableNames,faceData,streamData);
  }

  // parse which variables to average
  int nAvg = pp.countval("avgComps");
  Vector<std::string> avgComps;
  Vector<int> avgIdx(nAvg);
  if (nAvg>0) {
    pp.getarr("avgComps",avgComps);
    Print() << "Average components:" << std::endl;
    for (int iAvg=0; iAvg<nAvg; iAvg++) {
      avgIdx[iAvg] = getVarIdx(avgComps[iAvg],variableNames,nComps);
      Print() << "   " << avgIdx[iAvg] << ": " << avgComps[iAvg] << std::endl;
    }
  }

  // parse which variables to integrate
  int nInt = pp.countval("intComps");
  Vector<std::string> intComps;
  Vector<int> intIdx(nInt);
  if (nInt>0) {
    pp.getarr("intComps",intComps);
    Print() << "Integral components:" << std::endl;
    for (int iInt=0; iInt<nInt; iInt++) {
      intIdx[iInt] = getVarIdx(intComps[iInt],variableNames,nComps);
      Print() << "   " << intIdx[iInt] << ": " << intComps[iInt] << std::endl;
    }
  }

  int localDerOut = -1;
  // parse derived quantities (e.g. principal curvature zone)
  int nDerFlag = pp.countval("derComps");
  Vector<std::string> derCompsIn;
  Vector<std::string> derCompsOut;
  Vector<Vector<int>> derIdxIn(nDerFlag); //derIn
  Vector<Vector<int>> derIdxOut(nDerFlag); //derOut
  int nDerOut = 0;
  Real deltaT,rhoY,pkzFF;
  Real percOfMean(-1.);
  Real maxVolFac(-1.); // used to cap the integrals; a value of 3 may work well
  if (nDerFlag>0) pp.getarr("derComps",derCompsIn);
  for (int iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
    if (derCompsIn[iDerFlag]=="flameThickness") {
      Real reacTemp, prodTemp;
      pp.get("reacTemp",reacTemp);
      pp.get("prodTemp",prodTemp);
      deltaT=prodTemp-reacTemp;
      std::string tempGradVar="ModGradTemp";
      pp.query("tempGradVar",tempGradVar);
      derIdxIn[iDerFlag].push_back(getVarIdx(tempGradVar,variableNames,nComps)); //add gradT idx
      derIdxOut[iDerFlag].resize(1); //one out idx (thermalthickness)
      derIdxOut[iDerFlag][0] = (++localDerOut);
      derCompsOut.push_back("thermalThickness");
      nDerOut += 1;
    }
    if (derCompsIn[iDerFlag]=="flameSpeed") {
      pp.query("maxVolFac",maxVolFac);
      pp.query("percOfMean",percOfMean);
      pp.get("rhoY",rhoY);
      std::string FCRVar=fuelName+"_ConsumptionRate";
      pp.query("FCRVar",FCRVar);
      derIdxIn[iDerFlag].push_back(getVarIdx(FCRVar,variableNames,nComps)); //add FCR idx
      derIdxOut[iDerFlag].resize(2); //two out idx (flamespeed, reducedvol)
      for (int iDer = 0; iDer<2; iDer++) {
	derIdxOut[iDerFlag][iDer] = (++localDerOut);
      }
      derCompsOut.push_back("flameSpeed");
      derCompsOut.push_back("reducedVol");
      nDerOut += 2;
    }
    if (derCompsIn[iDerFlag]=="principalCurvatureZones") {
      std::string pkzMkVar = "MeanCurvature_prog_"+fuelName;
      std::string pkzGkVar = "GaussianCurvature_prog_"+fuelName;
      pp.query("pkzMkVar",pkzMkVar);
      pp.query("pkzGkVar",pkzGkVar);
      Real pkzLength;
      pp.get("pkzLength",pkzLength);
      pkzFF = 1.0 / (2*pkzLength);
      derIdxIn[iDerFlag].push_back(getVarIdx(pkzMkVar,variableNames,nComps)); //add mk idx
      derIdxIn[iDerFlag].push_back(getVarIdx(pkzGkVar,variableNames,nComps)); //add gk idx
      derIdxOut[iDerFlag].resize(3); //three out idx (zone,k1,k2)
      for (int iDer = 0; iDer<3; iDer++) {
	derIdxOut[iDerFlag][iDer] = (++localDerOut);
      }
      derCompsOut.push_back("k1");
      derCompsOut.push_back("k2");
      derCompsOut.push_back("zone");
      nDerOut += 3;      
    }
    if (derCompsIn[iDerFlag]=="reactionZoneThickness") {
      std::string reactionVar = "HeatRelease";
      pp.query("rztVar",reactionVar);
      derIdxIn[iDerFlag].push_back(getVarIdx(reactionVar,variableNames,nComps)); //add HRR idx
      derIdxOut[iDerFlag].resize(1); //one out idx (rzt)
      derIdxOut[iDerFlag][0] = (++localDerOut);
      derCompsOut.push_back("reactionZoneThickness");
      nDerOut += 1;
    }
  }
  AMREX_ALWAYS_ASSERT(derCompsOut.size() == nDerOut); //just a check on derive components out
  Print() << "Derived components: " << std::endl;
  for (int iDer = 0; iDer<nDerOut; iDer++) {
    Print() << derCompsOut[iDer] << std::endl;
  }
  
  // make space to hold output surface
  // for each element (i.e. triangle), we have three coordinates and one set of data
  // data will be written in triplicate, but no need to store all that
  // connectivity follows naturally from construction
  Vector<Real>         eltArea(nElts);
  Vector<Real>         eltVol(nElts);
  Vector<Real>         eltRatio(nElts);
  Vector<Array<dim3,AMREX_SPACEDIM>> surfLocs(nElts);
  Vector<Vector<Real>> surfAvg(nElts);
  Vector<Vector<Real>> surfInt(nElts);
  Vector<Vector<Real>> surfDer(nElts);
  Vector<int3> sIdx(nElts);
  // the surface is the mid point of the stream
  int surfPt = (nPtsOnStream-1)/2; // stream data location counts from zero    
  // set locations
  Print() << "Setting locations, making triangles, resizing arrays ..." << std::endl;
  
  // get the three stream indices from connectivity faceData
#ifdef _OPENMP
#pragma omp parallel for
#endif    
  for (int iElt=0; iElt<nElts; iElt++) {
    surfAvg[iElt].resize(nAvg);
    surfInt[iElt].resize(nInt);
    surfDer[iElt].resize(nDerOut);
    // define the spatial location of the two/three corners of the triangle
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) { 
      sIdx[iElt][iCorner] = faceData[iElt*AMREX_SPACEDIM+iCorner];
      for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	surfLocs[iElt][iCorner][d] = streamData[sIdx[iElt][iCorner]][nPtsOnStream*d+surfPt];
      }
    }
  }

  // stuff to dump PKZ streams
  int dumpPKZstreams = 0; pp.query("dumpPKZstreams", dumpPKZstreams);
  Vector<Array<Vector<Real>,6>> PKZstreams;
  Array<Real,6> areaZoneStreams;
  Array<Real,6> lsStreams;
  if(dumpPKZstreams) {
    PKZstreams.resize(nComps-AMREX_SPACEDIM+1);
    for (int i=0; i<nComps-AMREX_SPACEDIM+1; i++) {
      for (int j=0; j<6; j++) {
        PKZstreams[i][j].resize(nPtsOnStream);
        for (int k=0; k<nPtsOnStream; k++) {
          PKZstreams[i][j][k] = 0.0;
        }
      }
    }
    for (int i=0; i<6;i++) {
      areaZoneStreams[i] = 0;
      lsStreams[i] = 0.0;
    }
  }

  // evaluate area
  Print() << "Evaluating areas and volumes..." << std::endl;
  Real surfaceArea=0.;  
  Real totalVol=0.0;
  Real meanRatio = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: surfaceArea,totalVol) 
#endif
  for (int iElt=0; iElt<nElts; iElt++) {      
    // set up triangle ABC
    Array<dim3,AMREX_SPACEDIM> elt1,elt2;
    for (int d1=0; d1<AMREX_SPACEDIM; d1++) { // three components of location
      for (int d2=0; d2<AMREX_SPACEDIM; d2++) {
	elt1[d1][d2] = streamData[sIdx[iElt][d1]][nPtsOnStream*d2+surfPt];
      }
    }
    // find area
    eltArea[iElt] = elt_area(elt1); //(A,B,C)
    // keep a running total
    surfaceArea+=eltArea[iElt];
    
    eltVol[iElt]=0.;
    for (int iPt=1; iPt<nPtsOnStream; iPt++) {
      for (int d1=0; d1<AMREX_SPACEDIM; d1++) { // three components of location
	for (int d2=0; d2<AMREX_SPACEDIM; d2++) {
	  elt1[d1][d2] = streamData[sIdx[iElt][d1]][nPtsOnStream*d2+iPt-1];
	  elt2[d1][d2] = streamData[sIdx[iElt][d1]][nPtsOnStream*d2+iPt]; 
	}
      }
      Real vol = wedge_volume(elt1,elt2);
      eltVol[iElt] += vol;
      totalVol += vol;
    }
    eltRatio[iElt] = eltArea[iElt]/eltVol[iElt];
    meanRatio += eltArea[iElt]/eltVol[iElt];
  }
  meanRatio /= nElts;
  Print() << "   ... total surface area = " << surfaceArea << std::endl;
  Print() << "   ... total volume = " << totalVol << std::endl;

  // calculate averages
  //Print() << "Calculating averages ..." << std::endl;
  int iAvg, iInt, iDerFlag;
  Real filels = 0;
  Real filess = 0;
  Real filedelta = 0;
  Real fileEbar = 0;
  int getEbar = 0; pp.query("getEbar", getEbar);
  int strainIdx = -1;
  if(getEbar) {
    strainIdx = getVarIdx("StrainRate_prog_"+fuelName,avgComps,nAvg);
    amrex::Print() << "Getting Ebar, index: " << strainIdx << std::endl;
  }
  
  Real normAreaLS = surfaceArea; 
  Real normAreaSS = surfaceArea;
  Real normAreaEBAR = surfaceArea;
  Real normAreaDelta = surfaceArea;
  Print() << "Iterating over elements ..." << std::endl;
  int numFixedElts = 0;
  Real ds=-1.0;
  //calc streamLength from first stream of first element (should all be the same)?
  for (int iElt=0; iElt<nElts-1; iElt++) {
    Real s1 = 0.0;
    Real s2 = 0.0;
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      Real dsloc = (streamData[sIdx[iElt][0]][d*nPtsOnStream+1]-streamData[sIdx[iElt][0]][d*nPtsOnStream]);
      s1 += dsloc*dsloc;
      dsloc = (streamData[sIdx[iElt+1][0]][d*nPtsOnStream+1]-streamData[sIdx[iElt+1][0]][d*nPtsOnStream]);
      s2 += dsloc*dsloc;
    }
    s1 = std::sqrt(s1);
    s2 = std::sqrt(s2);
    if (s1 == s2) {
      ds = s1; //check later
      break;
    }
  }
  if (ds < 0) {
    Abort("ds broken");
  }
    
#ifdef _OPENMP
#pragma omp parallel for reduction(+: filess,filels,fileEbar,filedelta,numFixedElts)
#endif
  for (int iElt=0; iElt<nElts; iElt++) {
    //get thread-local stream index, data and element area
    int3 localSIdx = sIdx[iElt];
    Array<Vector<Real>,AMREX_SPACEDIM> localStreamData; 
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      localStreamData[d] = streamData[localSIdx[d]];
    }
    Real areaLoc = eltArea[iElt];
    Real volLoc = eltVol[iElt];
    Real ratioLoc = areaLoc/volLoc;//volLoc/areaLoc;

    // evaluate average of each component over the three points
    for (iAvg=0; iAvg<nAvg; iAvg++) {
      surfAvg[iElt][iAvg] = calcAvgVal(avgIdx[iAvg],nPtsOnStream,localStreamData);
    }
    if (getEbar) {
      if (std::isfinite(surfAvg[iElt][strainIdx])) {
	fileEbar += surfAvg[iElt][strainIdx]*areaLoc;
      } else {
	normAreaEBAR -= areaLoc;
      }
    }
    for (iInt=0; iInt<nInt; iInt++) {
      surfInt[iElt][iInt] = calcIntegral(intIdx[iInt],nPtsOnStream,localStreamData,areaLoc);
    }
    int outComp;
    for (iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
     if (derCompsIn[iDerFlag]=="flameThickness") {
       outComp = derIdxOut[iDerFlag][0];
       surfDer[iElt][outComp] = deltaT/calcMax(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData);
	// sum for mean
       if (std::isfinite(surfDer[iElt][outComp])) {
	 filels += surfDer[iElt][outComp] * areaLoc;
       } else {
	 normAreaLS -= areaLoc;
       }
     }
     if (derCompsIn[iDerFlag]=="reactionZoneThickness") {
       outComp=derIdxOut[iDerFlag][0];
       surfDer[iElt][outComp] = calcIntegral(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData,areaLoc)/calcMax(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData);
       if (std::isfinite(surfDer[iElt][outComp])) {
	 filedelta += surfDer[iElt][outComp]*areaLoc;
       } else {
	 normAreaDelta -= areaLoc;
       }
     }
     if (derCompsIn[iDerFlag]=="flameSpeed") {
       outComp = derIdxOut[iDerFlag][0];
       // cap the volume of element that can contribute to the integral
       if (maxVolFac>0.) {
	 surfDer[iElt][outComp] = calcCappedIntegral(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData,areaLoc,maxVolFac)/rhoY;
       } else {
	 surfDer[iElt][outComp] = calcIntegral(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData,areaLoc)/rhoY;
       }

       // sum for mean
       //Print() << ratioLoc/(ds*nPtsOnStream) << std::endl;
       if (std::isfinite(surfDer[iElt][outComp])) {
	 filess += surfDer[iElt][outComp] * areaLoc;
       } else {
	 normAreaSS -= areaLoc;
       }
     }
     if (derCompsIn[iDerFlag]=="principalCurvatureZones") {
       Real avgMk = calcAvgVal(derIdxIn[iDerFlag][0],nPtsOnStream,localStreamData);
       Real avgGk = calcAvgVal(derIdxIn[iDerFlag][1],nPtsOnStream,localStreamData);
       Real det = sqrt(fabs(avgMk*avgMk-avgGk));
       Real k1 = avgMk + det;
       Real k2 = avgMk - det;
       Real zone = 0;
       if   ( sqrt(k1*k1+k2*k2) <  pkzFF         )  zone = 1; // FF
       else {
	 if (                k2 >  half*k1       )  zone = 2; // LP
	 if (       fabs(k2/k1) <= half          )  zone = 3; // LE
	 if (        (-2*k1<k2) && (k2<-half*k1) )  zone = 4; // SP
	 if (       fabs(k1/k2) <= half          )  zone = 5; // TE
	 if (                k1 <  half*k2       )  zone = 6; // TP
       }
       if (zone == 0) {
	 Print() << "Could not find zone" << std::endl;
       }
       outComp = derIdxOut[iDerFlag][0];
       surfDer[iElt][outComp]=k1;
       outComp = derIdxOut[iDerFlag][1];
       surfDer[iElt][outComp]=k2;
       outComp = derIdxOut[iDerFlag][2];
       surfDer[iElt][outComp]=zone;
     }
    }

    // Thomas percOfMean correction
    for (iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
      if (derCompsIn[iDerFlag]=="flameSpeed") {
	if (percOfMean>0) {
	  for (int jDerFlag=0; jDerFlag<nDerFlag; jDerFlag++) {
	    if (derCompsIn[jDerFlag]=="principalCurvatureZones") {
	      Real zone = surfDer[iElt][derIdxOut[jDerFlag][2]];
	      if (zone < 1.5 || zone > 3.5) { //not leading point or edge
		if (ratioLoc < percOfMean*meanRatio) {
		  int nPtsOnReducedStream = (int)((ratioLoc/(percOfMean*meanRatio))*nPtsOnStream);
		  //Print() << nPtsOnReducedStream << "/" << nPtsOnStream << std::endl;
		  numFixedElts += 1;
		  if (nPtsOnReducedStream >= nPtsOnStream) {
		    Print() << nPtsOnReducedStream << std::endl;
		    Abort("Too many points on reduced stream");
		  }
		  if ( (nPtsOnStream - nPtsOnReducedStream) % 2 == 1) {
		    nPtsOnReducedStream -= 1;
		  }
		  outComp = derIdxOut[iDerFlag][0];
		  surfDer[iElt][outComp] = calcAdjustedIntegral(derIdxIn[iDerFlag][0],nPtsOnStream,nPtsOnReducedStream,localStreamData,areaLoc)/rhoY;
		  outComp = derIdxOut[iDerFlag][1];
		  surfDer[iElt][outComp] = 1;
		} else {
		  outComp = derIdxOut[iDerFlag][1];
		  surfDer[iElt][outComp] = 0;		   
		}
	      }
	    }
	  }
	}       
      }
    }
  
  }

  if (percOfMean>0) {
    Print() << "Number of elements length changed (area/volume ratio off): " << numFixedElts << "/" << nElts << std::endl;
    if (numFixedElts == nElts) {
      Print() << "ds = " << ds << std::endl;
      Abort("every element getting trimmed?");
    }
  }
  filels /= normAreaLS;
  filess /= normAreaSS;
  fileEbar /= normAreaEBAR;
  filedelta /= normAreaDelta;
  //array reducing stuff (JPDFs, conditionally averaged streamlines)
  //stuff for zone averaging

  Print() << "Sorting elements with out of bound streams" << std::endl;
  int tempIdx = getVarIdx("temp",variableNames,nComps);
  Vector<int> outOfBounds(nElts);
  int numOOB = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif    
  for (int iElt=0; iElt<nElts; iElt++) { 
    int3 localSIdx = sIdx[iElt];
    Array<Vector<Real>,AMREX_SPACEDIM> localStreamData; 
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      localStreamData[d] = streamData[localSIdx[d]];
    }
    dim3 endTemp, startTemp;
    //weird backwardness?
    for (int d=0;d<AMREX_SPACEDIM;d++) {
      endTemp[d] = localStreamData[d][nPtsOnStream*tempIdx]; 
      startTemp[d] = localStreamData[d][nPtsOnStream*tempIdx + nPtsOnStream - 1];
    }
    if (endTemp[0] <= startTemp[0] || endTemp[1] <= startTemp[1] || endTemp[2] <= startTemp[2]) {
      outOfBounds[iElt] = 1;
      numOOB += 1;
    } else {
      outOfBounds[iElt] = 0;
    }
  }

  Print() << "Number of elements with out of bound streams: " << numOOB << "/" << nElts << std::endl;
  
  int zoneComp = -1;
  int thermalThicknessComp = -1;
  int reactionZoneThicknessComp = -1;
  for (iDerFlag=0; iDerFlag<nDerFlag; iDerFlag++) {
    if (derCompsIn[iDerFlag]=="principalCurvatureZones") {
      zoneComp = derIdxOut[iDerFlag][2];
    }
    if (derCompsIn[iDerFlag]=="flameThickness") {
      thermalThicknessComp = derIdxOut[iDerFlag][0];
    }
    //if (derCompsIn[iDerFlag]=="reactionZoneThickness") {
    //  reactionZoneThicknessComp = derIdxOut[iDerFlag][0];
    //}
  }
  if (zoneComp < 0 && dumpPKZstreams) {
    Abort("can't find zone comp");
  }
  if (thermalThicknessComp < 0 && dumpPKZstreams) {
    Abort("can't find thermal thickness comp");
  }
  //if (reactionZoneThicknessComp < 0 && dumpPKZstream) {
    
  //}
#ifdef _OPENMP
  #pragma omp parallel
  {
#endif //declare local arrays here
    Vector<Array<Vector<Real>,6>> PKZstreams_local(PKZstreams);
    Array<Real,6> areaZoneStreams_local(areaZoneStreams);
    Array<Real,6> lsStreams_local(lsStreams);
#ifdef _OPENMP
#pragma omp for
#endif    
    for (int iElt=0; iElt<nElts; iElt++) { 
      int3 localSIdx = sIdx[iElt];
      Array<Vector<Real>,AMREX_SPACEDIM> localStreamData; 
      for (int d=0; d<AMREX_SPACEDIM; d++) {
	localStreamData[d] = streamData[localSIdx[d]];
      }
      
      Real areaLoc = eltArea[iElt];
      if (dumpPKZstreams && outOfBounds[iElt] == 0) {
	int zone = surfDer[iElt][zoneComp];
	areaZoneStreams_local[zone-1] += areaLoc;
	
	Real thermalThickness = surfDer[iElt][thermalThicknessComp];
	lsStreams_local[zone-1] += thermalThickness*areaLoc;
	dim3 surfPoint;
	for (int iComp=0; iComp<AMREX_SPACEDIM; iComp++) {
	  Real surfFaceVal = 0.0;
	  for (int d=0;d<AMREX_SPACEDIM;d++) {
	    surfFaceVal += localStreamData[d][nPtsOnStream*iComp + surfPt];
	  }
	  surfFaceVal /= static_cast<Real>(AMREX_SPACEDIM);
	  surfPoint[iComp] = surfFaceVal;
	}
	for (int iPt = 0; iPt<nPtsOnStream; iPt++) {	  
	  dim3 facePoint;   
	  for (int iComp=0; iComp<AMREX_SPACEDIM; iComp++) {
	    Real faceVal = 0.0;
	    for (int d=0;d<AMREX_SPACEDIM;d++) {
	      faceVal += localStreamData[d][nPtsOnStream*iComp + iPt];
	    }
	    faceVal /= static_cast<Real>(AMREX_SPACEDIM);
	    facePoint[iComp] = faceVal;	    
	  }
	  Real dist = 0.0;
	  for (int d=0;d<AMREX_SPACEDIM;d++) {
	    dist += (facePoint[d]-surfPoint[d])*(facePoint[d]-surfPoint[d]);
	  }
	  dist = std::sqrt(dist);
	  if (iPt > surfPt) {
	    dist *= -1;
	  }
	  
	  //Print() << zone << " " << iPt << std::endl;
	  PKZstreams_local[0][zone-1][iPt] += dist*areaLoc;
	}
	//probably swap these loop round but cba
	for (int iComp=AMREX_SPACEDIM; iComp<nComps; iComp++) {
	  for (int iPt = 0; iPt<nPtsOnStream; iPt++) {
	    Real faceVal = 0.0;
	    for (int d=0;d<AMREX_SPACEDIM;d++) {
	      faceVal += localStreamData[d][nPtsOnStream*iComp + iPt];
	    }
	    faceVal /= static_cast<Real>(AMREX_SPACEDIM);
	    PKZstreams_local[iComp-AMREX_SPACEDIM+1][zone-1][iPt] += faceVal*areaLoc;
	  }
	}
      }
    }
#ifdef _OPENMP
#pragma omp critical
  {
#endif
    if (dumpPKZstreams) {
      for (int z=0;z<6;z++) {
	areaZoneStreams[z] += areaZoneStreams_local[z];
	lsStreams[z] += lsStreams_local[z];
	for (int iComp=AMREX_SPACEDIM-1; iComp<nComps; iComp++) {
	  for (int iPt = 0; iPt<nPtsOnStream; iPt++) {
	    PKZstreams[iComp-AMREX_SPACEDIM+1][z][iPt] += PKZstreams_local[iComp-AMREX_SPACEDIM+1][z][iPt];
	  }
	}
      }
    }
#ifdef _OPENMP
  }
  }
#endif

  
  //dump characteristic values for this file
  std::string filename=infile+"/characteristics.dat";

  std::ofstream os(filename.c_str(),std::ios::out);
  
  os << filels << " " << filess << " " << fileEbar << " " << filedelta << std::endl;

  os.close();


  if (dumpPKZstreams) {
    std::string filename;
    Vector<std::string> zoneNames = {"FF","LP","LE","SP","TE","TP"};        

    for (int z = 0; z<6; z++) {
      filename = infile+"/"+zoneNames[z]+"_ls.dat";
      std::ofstream lsos(filename.c_str(),std::ios::out);
      if (areaZoneStreams[z] > 0) {
	lsStreams[z] /= areaZoneStreams[z];
      }
      lsos << lsStreams[z];
      lsos.close();

      filename = infile+"/"+zoneNames[z]+"_distfromseed.dat";
      std::ofstream dsos(filename.c_str(),std::ios::out);
      for (int iPt = nPtsOnStream-1; iPt>=0; iPt--) {
	if (areaZoneStreams[z] > 0) {
	  PKZstreams[0][z][iPt] /= areaZoneStreams[z];
	}
	dsos << PKZstreams[0][z][iPt] << " ";
      }
      dsos.close();
	
      for (int n = AMREX_SPACEDIM; n<nComps; n++) {
	filename = infile+"/"+zoneNames[z]+"_"+variableNames[n]+".dat";
	std::ofstream os(filename.c_str(),std::ios::out);
        for (int iPt = nPtsOnStream-1; iPt>=0; iPt--) {
	  if (areaZoneStreams[z] > 0) {
	    PKZstreams[n-AMREX_SPACEDIM+1][z][iPt] /= areaZoneStreams[z];
	  }
          os << PKZstreams[n-AMREX_SPACEDIM+1][z][iPt] << " ";
        }
        os.close();
      }
    }
  }

  if (writeSurface) {
    // write surface
  
    Print() << "Writing surface ..." << std::endl;
    writeSurfaceTecplot(infile,
			nElts, eltArea,  eltVol,  surfLocs,
			nAvg,  avgComps, surfAvg,
			nInt,  intComps, surfInt,
			nDerOut,  derCompsOut, surfDer);
    Print() << "   ... done" << std::endl;
  }
  // write basic
  if (writeBasic) {
    Print() << "Writing basic file ..." << std::endl;
    writeSurfaceBasic(infile,
		      nElts, eltArea,  eltVol,  surfLocs,
		      nAvg,  avgComps, surfAvg,
		      nInt,  intComps, surfInt,
		      nDerOut,  derCompsOut, surfDer);
    Print() << "   ... done" << std::endl;
  }

  return(0);
}

int getVarIdx(std::string varName, Vector<std::string>& variableNames, int nComps)
{  
    int varIdx=-1;
    for (int iComp=0; iComp<nComps; iComp++)
      if (variableNames[iComp]==varName) varIdx=iComp;
    if (varIdx==-1) {
      std::string msg=varName+" not in stream file";
      Abort(msg);
    }
    return varIdx;
}

Real calcAvgVal(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData)
{
  Real avgVal = 0.0;
  int surfPt = (nPtsOnStream-1)/2;
  // evaluate average of each component over the three points
  for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
    avgVal += streamData[iCorner][nPtsOnStream*compIdx+surfPt];
  }
  avgVal /= static_cast<Real>(AMREX_SPACEDIM);
  return avgVal;
}

Real calcMax(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData)
{
      Real maxVal = 0.0;
      for (int iPt=1; iPt<nPtsOnStream; iPt++) {
	Real avgVal = 0.0;
	for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	  avgVal += streamData[iCorner][nPtsOnStream*compIdx+iPt];
	}
	avgVal /= static_cast<Real>(AMREX_SPACEDIM);
	maxVal = max(maxVal,avgVal);
      }
      return maxVal;
}

//
// calculate the integral over the stream tube
//

Real calcIntegral(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData, Real eltArea)
{
  Real integral = 0.;
  
  // integrate
  Array<dim3,AMREX_SPACEDIM> elt1,elt2;
  dim3 val1, val2;
  for (int iPt=1; iPt<nPtsOnStream; iPt++) {
    for (int d1=0; d1<AMREX_SPACEDIM; d1++) { // three components of location
      for (int d2=0; d2<AMREX_SPACEDIM; d2++) {
	elt1[d1][d2] = streamData[d1][nPtsOnStream*d2+iPt-1];
	elt2[d1][d2] = streamData[d1][nPtsOnStream*d2+iPt];
      }
      val1[d1] = streamData[d1][nPtsOnStream*compIdx+iPt-1];
      val2[d1] = streamData[d1][nPtsOnStream*compIdx+iPt];
    }
    integral += wedge_volume_int(elt1,val1,elt2,val2);
  } 
  integral /= eltArea;
  return integral;
}

//
// calculate the integral, but cap the contribution by max volume
//

Real calcCappedIntegral(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData, Real eltArea, Real maxVolFac)
{
  Real integral = 0.;

  // integrate
  //dim3 A,B,C,D,E,F;
  Array<dim3,AMREX_SPACEDIM> elt1, elt2;
  dim3 val1,val2;
  // get reference volume at the surface
  // the surface is the mid point of the stream
  int surfPt = (nPtsOnStream-1)/2; // stream data location counts from zero
  // Now let's get the average volume of the two elements either side of the surface
  Real refVol(0.);
  for (int iPt=0; iPt<2; iPt++) {
    for (int d1=0; d1<AMREX_SPACEDIM; d1++) { // three components of location
      for (int d2=0; d2<AMREX_SPACEDIM; d2++) {
	elt1[d1][d2] = streamData[d1][nPtsOnStream*d2+surfPt+iPt-1];
	elt2[d1][d2] = streamData[d1][nPtsOnStream*d2+surfPt+iPt];
      }
    }
    refVol += 0.5*wedge_volume(elt1,elt2);
  }
  // set the maxVol to the factor passed in times this reference volume
  Real maxVol = maxVolFac*refVol;

  // now do the integral, capping by volFac
  for (int iPt=1; iPt<nPtsOnStream; iPt++) {
    for (int d1=0; d1<AMREX_SPACEDIM; d1++) { // three components of location
      for (int d2=0; d2<AMREX_SPACEDIM; d2++) {
	elt1[d1][d2] = streamData[d1][nPtsOnStream*d2+surfPt+iPt-1];
	elt2[d1][d2] = streamData[d1][nPtsOnStream*d2+surfPt+iPt];
      }
      val1[d1] = streamData[d1][nPtsOnStream*compIdx+iPt-1];
      val2[d1] = streamData[d1][nPtsOnStream*compIdx+iPt];           
    }

    // if volume of the element is bigger than the maximum reference volume,
    // then cap the contribution to the integral with volFac...
    // (should probably check that total consumption doesn't get missed)
    Real myVol = wedge_volume(elt1,elt2);
    Real volFac = min(myVol,maxVol)/(myVol+1.e-40);
    integral += volFac*wedge_volume_int(elt1,val1,elt2,val2);
  } 
  integral /= eltArea;
  return integral;
}

//
// calculate the integral on a reduced stream length
//


Real calcAdjustedIntegral(int compIdx, int nPtsOnStream, int nPtsOnReducedStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData, Real eltArea)
{
  Real integral = 0.;
  int diff = (nPtsOnStream - nPtsOnReducedStream)/2;
  // integrate
  Array<dim3,AMREX_SPACEDIM> elt1,elt2;
  dim3 val1,val2;
  for (int iPt=1; iPt<nPtsOnStream-diff; iPt++) {
    for (int d1=0; d1<AMREX_SPACEDIM; d1++) { // three components of location
      for (int d2=0; d2<AMREX_SPACEDIM; d2++) {
	elt1[d1][d2] = streamData[d1][nPtsOnStream*d2+iPt-1];
	elt2[d1][d2] = streamData[d1][nPtsOnStream*d2+iPt];
      }
      val1[d1] = streamData[d1][nPtsOnStream*compIdx+iPt-1];
      val2[d1] = streamData[d1][nPtsOnStream*compIdx+iPt];
    }
    integral += wedge_volume_int(elt1,val1,elt2,val2);
  } 
  integral /= eltArea;
  return integral;
}

//
// break up the variable list
//
std::vector<std::string> parseVarNames(std::istream& is)
{
    std::string line;
    std::getline(is,line);
    return amrex::Tokenize(line,std::string(" "));
}

//
// routine to read aja's binary stream files
//

void readStreamBin(std::string infile,
		   int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<int>& faceData,
		   Vector<Vector<Real>> &streamData)
{
  // Open header file
  std::string headerName = infile + "/Header";
  Print() << "Opening " << headerName << std::endl;
  std::ifstream ifs(headerName.c_str());
  std::istream* is = (infile=="-" ? (std::istream*)(&std::cin) : (std::istream*)(&ifs) );

  // read dummy header line
  std::string dummy;
  std::getline(ifs,dummy);

  // number of files to read
  int nFiles(-1);  
  ifs >> nFiles;
  Print() << "nFiles = " << nFiles << std::endl;

  // number of points on each stream
  ifs >> nStreams;
  Print() << "nStreams = " << nStreams << std::endl;

  // number of points on each stream
  ifs >> nPtsOnStream;
  Print() << "nPtsOnStream = " << nPtsOnStream << std::endl;

  // number of components
  ifs >> nComps;
  Print() << "nComps = " << nComps << std::endl;

  // next line
  std::getline(ifs,dummy);

  // read variable names
  variableNames.resize(nComps);
  variableNames = parseVarNames(*is);
  if (nComps!=variableNames.size())
    Abort("nComps != variableNames.size()");

  // connectivity data
  int fds;
  ifs >> fds;
  Print() << "faceDataSize = " << fds << std::endl;
  std::getline(ifs,dummy);
  faceData.resize(fds);
  ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());

  nElts = fds/3;
  Print() << "nElts = " << nElts << std::endl;

  // close header
  ifs.close();

  //
  // give the streams a home
  //
  streamData.resize(nStreams+1);
  for (int iStream=0; iStream<=nStreams; iStream++)
    streamData[iStream].resize(nPtsOnStream*nComps);
  
  // keep fread happy
  size_t read_size;

  //
  // now loop over binary files
  //
  for (int iFile=0; iFile<nFiles; iFile++) {
    //Print() << "Reading file " << iFile << std::endl;
    // read binary stream file
    std::string rootName = infile + "/str_";
    std::string fileName = Concatenate(rootName,iFile) + ".bin";
    FILE *file=fopen(fileName.c_str(),"r");

    int nFileStreams;
    read_size = fread(&(nFileStreams),sizeof(int),1,file);

    // loop over particle streams as written by pti (i.e. mangled order)
    for (int pindex=0; pindex<nFileStreams; pindex++) {
      // use the particle id to load data into right memory destination
      int iStream;
      read_size = fread(&(iStream),sizeof(int),1,file); // id
      
      for (int iComp=0; iComp<nComps; iComp++) {
	// by loading into [iStream] index, we're unmangling the parrallel particles
	int offset = iComp*nPtsOnStream;
	read_size = fread(&(streamData[iStream][offset]),sizeof(Real),nPtsOnStream,file);
      }      
    }
    
    fclose(file);

  } // iFiles

  Print() << "Finished reading stream data." << std::endl;
  
}

//
// write the stream data to matlab
//
void
writeStreamsMatlab(std::string infile,
		   int nStreams, int nPtsOnStream, int nComps,
		   std::vector<std::string>& variableNames,
		   Vector<Vector<Real>>& streamData)
{
  FILE *file;
  std::string filename;
  for (int iComp=0; iComp<nComps; iComp++) {
    filename = infile+"/"+variableNames[iComp]+".dat";
    file = fopen(filename.c_str(),"w");
    for (int iStream=1; iStream<=nStreams; iStream++) {
      for (int iPt=0; iPt<nPtsOnStream; iPt++) {
	fprintf(file,"%e ",streamData[iStream][nPtsOnStream*iComp+iPt]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  return;
}

//
// write the stream data (midpoint surface) straight to a tecplot file
//
void
writeSurfaceFromStreamTecplot(std::string infile,
			      int nStreams, int nElts, int nPtsOnStream, int nComps,
			      std::vector<std::string>& variableNames,
			      Vector<int>& faceData,
			      Vector<Vector<Real>>& streamData)
{
  std::string filename=infile+"_surfTec.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  std::string vars("VARIABLES =");
  for (int iComp=0; iComp<nComps; iComp++)
    vars += " " + variableNames[iComp];
  os << vars << std::endl;

  os << "ZONE T=\"streamBinTubeSurface\""
     << " N=" << nStreams
     << " E=" << nElts
     << " F=FEPOINT ET= TRIANGLE"
     << std::endl;

  int iPt=(nPtsOnStream-1)/2;
  for (int iStream=1; iStream<=nStreams; iStream++) {
    for (int iComp=0; iComp<nComps; iComp++) {
      os << streamData[iStream][nPtsOnStream*iComp+iPt] << " ";
    }
    os << std::endl;
  }

  for (int iElt=0; iElt<nElts; iElt++) {
    int offset=iElt*3;
    os << faceData[offset] << " "
       << faceData[offset+1] << " "
       << faceData[offset+2] << std::endl;
  }

  os.close();
  
  return;
}

//
// write all the surface quantities to a tecplot file
//
void
writeSurfaceTecplot(std::string infile,
		    int nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		    Vector<Array<dim3,AMREX_SPACEDIM>>& surfLocs,
		    int nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		    int nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		    int nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer)
{
  std::string filename=infile+"_binVolInt.dat";

  std::ofstream os(filename.c_str(),std::ios::out);
#if AMREX_SPACEDIM == 2
  std::string vars("VARIABLES = X Y area volume");
#else
  std::string vars("VARIABLES = X Y Z area volume");
#endif
  for (int iAvg=0; iAvg<nAvg; iAvg++)
    vars += " " + avgComps[iAvg] + "_avg";
  for (int iInt=0; iInt<nInt; iInt++)
    vars += " " + intComps[iInt] + "_volInt";
  for (int iDer=0; iDer<nDer; iDer++)
    vars += " " + derComps[iDer];
  os << vars << std::endl;

  os << "ZONE T=\"streamBinTubeSurface\""
     << " N=" << nElts*AMREX_SPACEDIM
     << " E=" << nElts
     << " F=FEPOINT ET=TRIANGLE"
     << std::endl;

  // write averages
  os << std::setprecision(12);
  for (int iElt=0; iElt<nElts; iElt++) {
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
      // coordinate
      for (int d=0; d<AMREX_SPACEDIM; d++)
	os << surfLocs[iElt][iCorner][d] << " ";
      // area
      os << eltArea[iElt] << " ";
      // volume
      os << eltVol[iElt] << " ";
      // averages
      for (int iAvg=0; iAvg<nAvg; iAvg++)
	os << surfAvg[iElt][iAvg] << " ";
      // integrals
      for (int iInt=0; iInt<nInt; iInt++)
	os << surfInt[iElt][iInt] << " ";
      // derived
      for (int iDer=0; iDer<nDer; iDer++)
	os << surfDer[iElt][iDer] << " ";
      os << std::endl;
    }
  }

  // write connectivity
  for (int iElt=1; iElt<3*nElts;)
    os << iElt++ << " " << iElt++ << " " << iElt++ << std::endl;

  os.close();
}

//
// write all the surface quantities to a tecplot file
//
void
writeSurfaceBasic(std::string infile,
		  int nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		  Vector<Array<dim3,AMREX_SPACEDIM>>& surfLocs,
		  int nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		  int nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		  int nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer)
{
  std::string filename=infile+"_binVolInt_basic.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  // write averages
  os << std::setprecision(12);
  for (int iElt=0; iElt<nElts; iElt++) {
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
      // coordinate
      for (int d=0; d<AMREX_SPACEDIM; d++)
	os << surfLocs[iElt][iCorner][d] << " ";
      // area
      os << eltArea[iElt] << " ";
      // volume
      os << eltVol[iElt] << " ";
      // averages
      for (int iAvg=0; iAvg<nAvg; iAvg++)
	os << surfAvg[iElt][iAvg] << " ";
      // integrals
      for (int iInt=0; iInt<nInt; iInt++)
	os << surfInt[iElt][iInt] << " ";
      // derived
      for (int iDer=0; iDer<nDer; iDer++)
	os << surfDer[iElt][iDer] << " ";
      os << std::endl;
    }
  }

  os.close();
}

#if AMREX_SPACEDIM == 2
// -----------------------------------------------------------------------------
// Helper: area of triangle (2D)
// -----------------------------------------------------------------------------
Real triArea(const dim3& A, const dim3& B, const dim3& C)
{
  return half * std::abs((B[0]-A[0])*(C[1]-A[1]) - (C[0]-A[0])*(B[1]-A[1]));
}
#else
// -----------------------------------------------------------------------------
// Helper: volume of tetrahedron (3D)
// -----------------------------------------------------------------------------
Real tetVol(const dim3& A, const dim3& B,
	    const dim3& C, const dim3& D)
{
  dim3 V1, V2, V3, cross;
  for (int i=0;i<3;++i) {
    V1[i] = B[i]-A[i];
    V2[i] = C[i]-A[i];
    V3[i] = D[i]-A[i];
  }
  cross[0] = V2[1]*V3[2] - V2[2]*V3[1];
  cross[1] = V2[2]*V3[0] - V2[0]*V3[2];
  cross[2] = V2[0]*V3[1] - V2[1]*V3[0];
  Real vol = (V1[0]*cross[0] + V1[1]*cross[1] + V1[2]*cross[2])/6.0;
  return std::abs(vol);
}
#endif


Real wedge_volume(const Array<dim3,AMREX_SPACEDIM>& elt1,
                  const Array<dim3,AMREX_SPACEDIM>& elt2)
{
#if AMREX_SPACEDIM == 2
  // 2D: quadrilateral / parallelogram
  // elt1 = [A,B], el 2 = [C,D]  (N==2)
  const dim3& A = elt1[0];
  const dim3& B = elt1[1];
  const dim3& C = elt2[0];
  const dim3& D = elt2[1];	
  return triArea(A,B,C) + triArea(A,C,D);
#else
  // 3D: triangular prism / wedge
  // elt1 = [A,B,C], elt2 = [D,E,F]
  const dim3& A = elt1[0];
  const dim3& B = elt1[1];
  const dim3& C = elt1[2];
  const dim3& D = elt2[0];
  const dim3& E = elt2[1];
  const dim3& F = elt2[2];
  
  return tetVol(A,B,C,E) + tetVol(A,D,E,F) + tetVol(A,C,E,F);
#endif
}

Real elt_area(const Array<dim3,AMREX_SPACEDIM>& elt) {
#if AMREX_SPACEDIM == 2
  // Line segment length
  Real sum = 0;
  for (size_t i = 0; i < AMREX_SPACEDIM; ++i) {
    Real d = elt[1][i] - elt[0][i];
    sum += d*d;
  }
  return std::sqrt(sum);
#else
  // Triangle area
  dim3 R1, R2, R3;

  for (int i=0; i<AMREX_SPACEDIM; ++i) {
    R1[i] = elt[1][i] - elt[0][i]; // B-A
    R2[i] = elt[2][i] - elt[0][i]; // C-A
  }
  
  // Cross product
  R3[0] = R1[1]*R2[2] - R2[1]*R1[2];
  R3[1] = R1[2]*R2[0] - R2[2]*R1[0];
  R3[2] = R1[0]*R2[1] - R2[0]*R1[1];
  
  Real result = 0;
  for (int i=0; i<3; ++i) {
    result += R3[i]*R3[i];
  }
  return 0.5 * std::sqrt(result);
#endif
}



// -----------------------------------------------------------------------------
// Generic wedge_volume_int for 2D quadrilateral or 3D wedge
//   - DIM = 2 → quadrilateral (A,B in first row; C,D in second row)
//   - DIM = 3 → wedge (A,B,C in first row; D,E,F in second row)
// -----------------------------------------------------------------------------
Real wedge_volume_int(const Array<dim3,AMREX_SPACEDIM>& elt1,
                      const dim3& val1,
                      const Array<dim3,AMREX_SPACEDIM>& elt2,
                      const dim3& val2)
{
#if AMREX_SPACEDIM == 2
  const dim3& A = elt1[0];
  const dim3& B = elt1[1];
  const dim3& C = elt2[0];
  const dim3& D = elt2[1];
  
  Real vA = val1[0], vB = val1[1];
  Real vC = val2[0], vD = val2[1];
  
  // Sub-areas like sub-tetrahedra in 3D
  Real area_ABC = triArea(A,B,C);
  Real area_ACD = triArea(A,C,D);
  Real area_ABD = triArea(A,B,D);
  Real area_BCD = triArea(B,C,D);
  
  // Integrals over sub-triangles
  Real int_1 = (vA+vB+vC) * area_ABC / 3.0;
  Real int_2 = (vA+vC+vD) * area_ACD / 3.0;
  Real int_3 = (vA+vB+vD) * area_ABD / 3.0;
  Real int_4 = (vB+vC+vD) * area_BCD / 3.0;
  
  // Average contributions for higher-order accuracy
  return quarter * (int_1 + int_2 + int_3 + int_4);
#else
  // 3D wedge: A,B,C bottom; D,E,F top
  const dim3& A = elt1[0]; const Real vA = val1[0];
  const dim3& B = elt1[1]; const Real vB = val1[1];
  const dim3& C = elt1[2]; const Real vC = val1[2];
  const dim3& D = elt2[0]; const Real vD = val2[0];
  const dim3& E = elt2[1]; const Real vE = val2[1];
  const dim3& F = elt2[2]; const Real vF = val2[2];
  
  // replicate old trusted method
  const Real vol_EABC = tetVol(A,B,C,E);
  const Real vol_ADEF = tetVol(A,D,E,F);
  const Real vol_ACEF = tetVol(C,E,F,A);
  const Real vol_DABC = tetVol(A,B,C,D);
  const Real vol_FABC = tetVol(A,B,C,F);
  const Real vol_BDEF = tetVol(B,D,E,F);
  const Real vol_CDEF = tetVol(C,D,E,F);
  const Real vol_ACED = tetVol(C,E,D,A);
  const Real vol_BCDF = tetVol(B,C,D,F);
  const Real vol_BCDE = tetVol(B,C,D,E);
  const Real vol_ABDF = tetVol(B,D,F,A);
  const Real vol_ABEF = tetVol(B,E,F,A);
  
  const Real int_1 = ((vD+vA+vB+vC)*vol_DABC +
		      (vB+vD+vE+vF)*vol_BDEF +
		      (vB+vC+vD+vF)*vol_BCDF)*quarter;
  
  const Real int_2 = ((vD+vA+vB+vC)*vol_DABC +
		      (vC+vD+vE+vF)*vol_CDEF +
		      (vB+vC+vD+vE)*vol_BCDE)*quarter;
  
  const Real int_3 = ((vE+vA+vB+vC)*vol_EABC +
		      (vA+vD+vE+vF)*vol_ADEF +
		      (vA+vC+vE+vF)*vol_ACEF)*quarter;
  
  const Real int_4 = ((vE+vA+vB+vC)*vol_EABC +
		      (vC+vD+vE+vF)*vol_CDEF +
		      (vA+vC+vE+vD)*vol_ACED)*quarter;
  
  const Real int_5 = ((vF+vA+vB+vC)*vol_FABC +
		      (vA+vD+vE+vF)*vol_ADEF +
		      (vA+vB+vE+vF)*vol_ABEF)*quarter;
  
  const Real int_6 = ((vF+vA+vB+vC)*vol_FABC +
		      (vB+vD+vE+vF)*vol_BDEF +
		      (vA+vB+vD+vF)*vol_ABDF)*quarter;
  
  return (int_1+int_2+int_3+int_4+int_5+int_6)*sixth;
#endif
}

