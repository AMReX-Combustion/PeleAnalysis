#include <string>
#include <iostream>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

#include <streamBinTubeStats.H>

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  if (ParallelDescriptor::NProcs()>1)
    Abort("Code is not yet parallel safe");
  
  ParmParse pp;

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
  
  // declare size and data holders
  int nStreams, nElts, nPtsOnStream, nComps;
  std::vector<std::string> variableNames;
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
      avgIdx[iAvg]=-1;
      for (int iComp=0; iComp<nComps; iComp++)
	if (variableNames[iComp]==avgComps[iAvg]) avgIdx[iAvg]=iComp;
      if (avgIdx[iAvg]==-1) {
	std::string msg="avgComp (" + avgComps[iAvg] + ") not in stream file";
	Abort(msg);
      }
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
      intIdx[iInt]=-1;
      for (int iComp=0; iComp<nComps; iComp++)
	if (variableNames[iComp]==intComps[iInt]) intIdx[iInt]=iComp;
      if (intIdx[iInt]==-1) {
	std::string msg="intComp (" + intComps[iInt] + ") not in stream file";
	Abort(msg);
      }
      Print() << "   " << intIdx[iInt] << ": " << intComps[iInt] << std::endl;
    }
  }

  // parse derived quantities (e.g. principal curvature zone)
  int nDer = pp.countval("derComps");
  Vector<std::string> derComps;
  if (nDer>0) pp.getarr("derComps",derComps);
  
  // make space to hold output surface
  // for each element (i.e. triangle), we have three coordinates and one set of data
  // data will be written in triplicate, but no need to store all that
  // connectivity follows naturally from construction
  Vector<Real>         eltArea(nElts);
  Vector<Real>         eltVol(nElts);
  Vector<Vector<dim3>> surfLocs(nElts);
  Vector<Vector<Real>> surfAvg(nElts);
  Vector<Vector<Real>> surfInt(nElts);
  Vector<Vector<Real>> surfDer(nElts);
  for (int iElt=0; iElt<nElts; iElt++) {
    surfLocs[iElt].resize(AMREX_SPACEDIM);
    surfAvg[iElt].resize(nAvg);
    surfInt[iElt].resize(nInt);
    surfDer[iElt].resize(nDer);
  }
 
  // get the three stream indices from connectivity faceData
  Print() << "Making triangles from connectivity data ..." << std::endl;
  Vector<int3> sIdx(nElts);
  for (int iElt=0; iElt<nElts; iElt++)
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
      sIdx[iElt][iCorner] = faceData[iElt*AMREX_SPACEDIM+iCorner];
  
  // the surface is the mid point of the stream
  int surfPt = (nPtsOnStream-1)/2; // stream data location counts from zero

  // set locations
  Print() << "Setting locations ..." << std::endl;
  for (int iElt=0; iElt<nElts; iElt++) {
    // define the spatial location of the three corners of the triangle
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) { // three corners (streams)
      int iStream=faceData[iElt*AMREX_SPACEDIM+iCorner];
      for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	surfLocs[iElt][iCorner][d] = streamData[iStream][nPtsOnStream*d+surfPt];
      }
    }
  }

  // evaluate area
  Print() << "Evaluating areas ..." << std::endl;
  Real surfaceArea=0.;
  for (int iElt=0; iElt<nElts; iElt++) {
    int3 iStream; // stream indices that make up this element
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
      iStream[iCorner]=sIdx[iElt][iCorner];
    // set up triangle ABC
    dim3 A, B, C;
    for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
      A[d] = streamData[iStream[0]][nPtsOnStream*d+surfPt];
      B[d] = streamData[iStream[1]][nPtsOnStream*d+surfPt];
      C[d] = streamData[iStream[2]][nPtsOnStream*d+surfPt];
    }
    // find area
    eltArea[iElt] = wedge_area(A,B,C);
    // keep a running total
    surfaceArea+=eltArea[iElt];
  }
  Print() << "   ... total surface area = " << surfaceArea << std::endl;

  // evaluate volume
  Print() << "Evaluating volumes ..." << std::endl;
  for (int iElt=0; iElt<nElts; iElt++) {
    int3 iStream; // stream indices that make up this element
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
      iStream[iCorner]=sIdx[iElt][iCorner];
    eltVol[iElt]=0.;
    for (int iPt=1; iPt<nPtsOnStream; iPt++) {
      dim3 A, B, C, D, E, F;
      for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	A[d] = streamData[iStream[0]][nPtsOnStream*d+iPt-1];
	B[d] = streamData[iStream[1]][nPtsOnStream*d+iPt-1];
	C[d] = streamData[iStream[2]][nPtsOnStream*d+iPt-1];
	D[d] = streamData[iStream[0]][nPtsOnStream*d+iPt];
	E[d] = streamData[iStream[1]][nPtsOnStream*d+iPt];
	F[d] = streamData[iStream[2]][nPtsOnStream*d+iPt];
      }
      eltVol[iElt] += wedge_volume(A,B,C,D,E,F);
    }
  }

  // calculate averages
  Print() << "Calculating averages ..." << std::endl;
  for (int iElt=0; iElt<nElts; iElt++) {
    // evaluate average of each component over the three points
    for (int iAvg=0; iAvg<nAvg; iAvg++) {
      int iComp=avgIdx[iAvg]; // component index that we're averaging
      Real avgVal(0.);
      for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	int iStream=sIdx[iElt][iCorner];
	avgVal += streamData[iStream][nPtsOnStream*iComp+surfPt];
      }
      surfAvg[iElt][iAvg] = avgVal*third;
    }
  }
  
  // calculate integrals
  Print() << "Calculating integrals ..." << std::endl;
  for (int iElt=0; iElt<nElts; iElt++) {
    int3 iStream; // stream indices that make up this element
    for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) 
      iStream[iCorner]=sIdx[iElt][iCorner];
    // set integral to zero
    for (int iInt=0; iInt<nInt; iInt++) 
      surfInt[iElt][iInt] = 0.;
    // integrate
    for (int iPt=1; iPt<nPtsOnStream; iPt++) {
      dim3 A, B, C, D, E, F;
      for (int d=0; d<AMREX_SPACEDIM; d++) { // three components of location
	A[d] = streamData[iStream[0]][nPtsOnStream*d+iPt-1];
	B[d] = streamData[iStream[1]][nPtsOnStream*d+iPt-1];
	C[d] = streamData[iStream[2]][nPtsOnStream*d+iPt-1];
	D[d] = streamData[iStream[0]][nPtsOnStream*d+iPt];
	E[d] = streamData[iStream[1]][nPtsOnStream*d+iPt];
	F[d] = streamData[iStream[2]][nPtsOnStream*d+iPt];
      }
      
      for (int iInt=0; iInt<nInt; iInt++) {
	int iComp = intIdx[iInt];
	Real vA = streamData[iStream[0]][nPtsOnStream*iComp+iPt-1];
	Real vB = streamData[iStream[1]][nPtsOnStream*iComp+iPt-1];
	Real vC = streamData[iStream[2]][nPtsOnStream*iComp+iPt-1];
	Real vD = streamData[iStream[0]][nPtsOnStream*iComp+iPt];
	Real vE = streamData[iStream[1]][nPtsOnStream*iComp+iPt];
	Real vF = streamData[iStream[2]][nPtsOnStream*iComp+iPt];
	surfInt[iElt][iInt] += wedge_volume_int(A,vA,B,vB,C,vC,D,vD,E,vE,F,vF);
      }
    }
    for (int iInt=0; iInt<nInt; iInt++) 
      surfInt[iElt][iInt] /= eltArea[iElt];
  }
  
  // calculate derived
  Print() << "Calculating derived quanities ..." << std::endl;
  for (int iDer=0; iDer<nDer; iDer++) {

    //
    // local flame thickness
    //
    if (derComps[iDer]=="flameThickness") {
      Print() << "   Evaluating local flame thermal thickness..." << std::endl;
      
      // need reactant and product temperatures
      Real reacTemp, prodTemp;
      pp.get("reacTemp",reacTemp);
      pp.get("prodTemp",prodTemp);
      Real deltaT=prodTemp-reacTemp;
      
      // also need to know which variable to use as temperature gradient
      std::string tempGradVar;
      pp.get("tempGradVar",tempGradVar);
      int tempGradComp=-1;
      for (int iComp=0; iComp<nComps; iComp++)
	if (variableNames[iComp]==tempGradVar) tempGradComp=iComp;
      if (tempGradComp==-1) {
	std::string msg="tempGradComp (" + tempGradVar + ") not in stream file";
	Abort(msg);
      }
      Print() << "      " << tempGradComp << ": " << tempGradVar << std::endl
	      << "      " << "deltaT = " << prodTemp << " - " << reacTemp
	      << " = " << deltaT << std::endl;

      // mean local thermal thickness
      Real ls(0.);
      for (int iElt=0; iElt<nElts; iElt++) {
	Real maxModGradT(0.);
	for (int iPt=1; iPt<nPtsOnStream; iPt++) {
	  Real avgModGradT(0.);
	  for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	    int iStream=sIdx[iElt][iCorner];
	    avgModGradT += streamData[iStream][nPtsOnStream*tempGradComp+iPt];
	  }
	  avgModGradT *= third;
	  maxModGradT = max(maxModGradT,avgModGradT);
	}
	// local thermal thicnkess
	surfDer[iElt][iDer] = deltaT/maxModGradT;
	// sum for mean
	ls += surfDer[iElt][iDer] * eltArea[iElt];
      }
      ls /= surfaceArea;
      Print () << "      mean local thermal thickness = " << ls << std::endl;
    }

    //
    // principal curvature zone
    //
    if (derComps[iDer]=="principalCurvatureZone") {
      Print() << "   Evaluating principal curvature zone..." << std::endl;
      // need a length scale to define "flat"
      Real pkzLength;
      pp.get("pkzLength",pkzLength);
      Real pkzFF = half/pkzLength;

      // also need to know which variable to use as mean and gaussian curvature
      std::string pkzMkVar, pkzGkVar;
      pp.get("pkzMkVar",pkzMkVar);
      pp.get("pkzGkVar",pkzGkVar);
      int pkzMkComp=-1;
      int pkzGkComp=-1;
      for (int iComp=0; iComp<nComps; iComp++) {
	if (variableNames[iComp]==pkzMkVar) pkzMkComp=iComp;
	if (variableNames[iComp]==pkzGkVar) pkzGkComp=iComp;
      }
      if (pkzMkComp==-1) {
	std::string msg="pkzMkComp (" + pkzMkVar + ") not in stream file";
	Abort(msg);
      }
      if (pkzGkComp==-1) {
	std::string msg="pkzGkComp (" + pkzGkVar + ") not in stream file";
	Abort(msg);
      }
      Print() << "      " << pkzMkComp << ": " << pkzMkVar << std::endl
	      << "      " << pkzGkComp << ": " << pkzGkVar << std::endl
	      << "      " << "pkzLength = " << pkzLength << std::endl;

      for (int iElt=0; iElt<nElts; iElt++) {
	Real avgMk(0.);
	Real avgGk(0.);
	for (int iCorner=0; iCorner<AMREX_SPACEDIM; iCorner++) {
	  int iStream=sIdx[iElt][iCorner];
	  avgMk += streamData[iStream][nPtsOnStream*pkzMkComp+surfPt];
	  avgGk += streamData[iStream][nPtsOnStream*pkzGkComp+surfPt];
	}
	avgMk *= third;
	avgGk *= third;
	Real det = sqrt(fabs(avgMk*avgMk-avgGk));
	Real k1 = avgMk + det;
	Real k2 = avgMk - det;
	Real zone(0);
	if   ( sqrt(k1*k1+k2*k2) <  pkzFF         )  zone = 1; // FF
	else {
	  if (                k2 >  half*k1       )  zone = 2; // LP
	  if (       fabs(k2/k1) <= half          )  zone = 3; // LE
	  if (        (-2*k1<k2) && (k2<-half*k1) )  zone = 4; // SP
	  if (       fabs(k1/k2) <= half          )  zone = 5; // TE
	  if (                k1 <  half*k2       )  zone = 6; // TP
	}
	surfDer[iElt][iDer] = zone;
      }
    } // pkz
    
  }
  
  // write surface
  Print() << "Writing surface ..." << std::endl;
  writeSurfaceTecplot(infile,
		      nElts, eltArea,  eltVol,  surfLocs,
		      nAvg,  avgComps, surfAvg,
		      nInt,  intComps, surfInt,
		      nDer,  derComps, surfDer);
  Print() << "   ... done" << std::endl;

  // write basic
  if (writeBasic) {
    Print() << "Writing basic file ..." << std::endl;
    writeSurfaceBasic(infile,
		      nElts, eltArea,  eltVol,  surfLocs,
		      nAvg,  avgComps, surfAvg,
		      nInt,  intComps, surfInt,
		      nDer,  derComps, surfDer);
    Print() << "   ... done" << std::endl;
  }

  return(0);
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
		   int& nStreams, int& nPtsOnStream, int& nComps,
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
			      int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
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
		    int& nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		    Vector<Vector<dim3>>& surfLocs,
		    int& nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		    int& nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		    int& nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer)
{
  std::string filename=infile+"_binVolInt.dat";

  std::ofstream os(filename.c_str(),std::ios::out);

  std::string vars("VARIABLES = X Y Z area volume");
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
		  int& nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		  Vector<Vector<dim3>>& surfLocs,
		  int& nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		  int& nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		  int& nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer)
{
  std::string filename=infile+"_binVolInt.dat";

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

//
// area of the triangle
//
Real
wedge_area(const dim3& A,
	   const dim3& B,
	   const dim3& C)
{
  Real result;  
  Real R1[3], R2[3], R3[3];
  
  for (int i=0; i<3; ++i) {
    R1[i] = B[i] - A[i];
    R2[i] = C[i] - A[i];
  }
  
  R3[0] = R1[1]*R2[2] - R2[1]*R1[2];
  R3[1] = R1[2]*R2[0] - R2[2]*R1[0];
  R3[2] = R1[0]*R2[1] - R2[0]*R1[1];

  result=0;
  for (int i=0; i<3; ++i) {
    result += R3[i]*R3[i];
  }
  result = half*std::sqrt(result);

  return result;
}

//
// volume of a tetrahedron
//
Real
tetVol(const dim3& A, const dim3& B, const dim3& C, const dim3& D)
{
    Real R1[3], R2[3], R3[3], R4[3];

    for (int i=0; i<3; ++i)
    {
        R1[i] = D[i] - A[i];
        R2[i] = B[i] - A[i];
        R3[i] = C[i] - A[i];
    }

    for (int i=0; i<3; ++i)
    {
        R4[0] = R2[1]*R3[2] - R3[1]*R2[2];
        R4[1] = R2[2]*R3[0] - R3[2]*R2[0];
        R4[2] = R2[0]*R3[1] - R3[0]*R2[1];
    }
    
    Real res=0;
    for (int i=0; i<3; ++i)
        res += R1[i]*R4[i];

    return std::abs(res*sixth);
}

//
// volume of the irregular triangular prism
//
Real
wedge_volume(const dim3& A, const dim3& B, const dim3& C,
	     const dim3& D, const dim3& E, const dim3& F)
{
  Real result;
  
  const Real vol_EABC = tetVol(A,B,C,E);
  const Real vol_ADEF = tetVol(A,D,E,F);
  const Real vol_ACEF = tetVol(C,E,F,A);
  
  result = (vol_EABC + vol_ADEF + vol_ACEF);
  
  return result;
}

//
// integral over the irregular triangular prism
//
Real
wedge_volume_int(const dim3& A, const Real vA,
		 const dim3& B, const Real vB,
		 const dim3& C, const Real vC,
		 const dim3& D, const Real vD,
		 const dim3& E, const Real vE,
		 const dim3& F, const Real vF)
{
  Real result;

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
  
  const Real int_1 = ( (vD+vA+vB+vC)*vol_DABC + 
		       (vB+vD+vE+vF)*vol_BDEF + 
		       (vB+vC+vD+vF)*vol_BCDF )*quarter;
  
  const Real int_2 = ( (vD+vA+vB+vC)*vol_DABC + 
		       (vC+vD+vE+vF)*vol_CDEF + 
		       (vB+vC+vD+vE)*vol_BCDE )*quarter;
  
  const Real int_3 = ( (vE+vA+vB+vC)*vol_EABC + 
		       (vA+vD+vE+vF)*vol_ADEF + 
		       (vA+vC+vE+vF)*vol_ACEF )*quarter;
  
  const Real int_4 = ( (vE+vA+vB+vC)*vol_EABC + 
		       (vC+vD+vE+vF)*vol_CDEF + 
		       (vA+vC+vE+vD)*vol_ACED )*quarter;
  
  const Real int_5 = ( (vF+vA+vB+vC)*vol_FABC + 
		       (vA+vD+vE+vF)*vol_ADEF + 
		       (vA+vB+vE+vF)*vol_ABEF )*quarter;
  
  const Real int_6 = ( (vF+vA+vB+vC)*vol_FABC + 
		       (vB+vD+vE+vF)*vol_BDEF + 
		       (vA+vB+vD+vF)*vol_ABDF )*quarter;
  
  result = (int_1 + int_2 + int_3 + int_4 + int_5 + int_6)*sixth;
  
  return result;
}

