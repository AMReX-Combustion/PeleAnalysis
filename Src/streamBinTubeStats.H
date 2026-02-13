typedef Array<int, AMREX_SPACEDIM> int3;
typedef Array<Real,AMREX_SPACEDIM> dim3;

const Real half    = 1./2.;
const Real third   = 1./3.;
const Real quarter = 1./4.;
const Real sixth   = 1./6.;

std::vector<std::string> parseVarNames(std::istream& is);

int getVarIdx(std::string varName, Vector<std::string>& variableNames, int nComps);

Real calcAvgVal(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData);

Real calcMax(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData);

Real calcIntegral(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData, Real eltArea, const Vector<int>& is_per_dim, const Array<Real,AMREX_SPACEDIM>& domain_size);

Real calcCappedIntegral(int compIdx, int nPtsOnStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData, Real eltArea, Real maxVolFac, const Vector<int>& is_per_dim, const Array<Real,AMREX_SPACEDIM>& domain_size);

Real calcAdjustedIntegral(int compIdx, int nPtsOnStream, int nPtsOnReducedStream, Array<Vector<Real>,AMREX_SPACEDIM>& streamData, Real eltArea, const Vector<int>& is_per_dim, const Array<Real,AMREX_SPACEDIM>& domain_size);

//int streamReduce(int nPtsOnStream, Vector<Vector<Real>>& streamData, Real eltArea, Real eltVol);

void readStreamBin(std::string infile,
		   int& nStreams, int& nElts, int& nPtsOnStream, int& nComps,
		   std::vector<std::string>& variableNames,
		   Vector<int>& faceData,
		   Vector<Vector<Real>> &streamData);

void
writeStreamsMatlab(std::string infile,
		   int nStreams, int nPtsOnStream, int nComps,
		   std::vector<std::string>& variableNames,
		   Vector<Vector<Real>>& streamData);
  
void
writeSurfaceFromStreamTecplot(std::string infile,
			      int nStreams, int nElts, int nPtsOnStream, int nComps,
			      std::vector<std::string>& variableNames,
			      Vector<int>& faceData,
			      Vector<Vector<Real>>& streamData);
  
void
writeSurfaceTecplot(std::string infile,
		    int nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		    Vector<Array<dim3,AMREX_SPACEDIM>>& surfLocs,
		    int nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		    int nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		    int nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer);

void
writeSurfaceBasic(std::string infile,
		  int nElts, Vector<Real>&        eltArea,  Vector<Real>&        eltVol,
		  Vector<Array<dim3,AMREX_SPACEDIM>>& surfLocs,
		  int nAvg,  Vector<std::string>& avgComps, Vector<Vector<Real>>& surfAvg,
		  int nInt,  Vector<std::string>& intComps, Vector<Vector<Real>>& surfInt,
		  int nDer,  Vector<std::string>& derComps, Vector<Vector<Real>>& surfDer);

Real
elt_area(const Array<dim3,AMREX_SPACEDIM>& elt, const Vector<int>& is_per_dim, const Array<Real,AMREX_SPACEDIM>& domain_size);

#if AMREX_SPACEDIM == 2
// -----------------------------------------------------------------------------
// Helper: area of triangle (2D)
// -----------------------------------------------------------------------------
Real triArea(const dim3& A, const dim3& B, const dim3& C, 
	         const Vector<int>& is_per_dim, 
			 const Array<Real,AMREX_SPACEDIM>& domain_size);
#elif AMREX_SPACEDIM == 3
// -----------------------------------------------------------------------------
// Helper: volume of tetrahedron (3D)
// -----------------------------------------------------------------------------
Real tetVol(const dim3& A, const dim3& B, const dim3& C, const dim3& D, 
	        const Vector<int>& is_per_dim, 
			const Array<Real,AMREX_SPACEDIM>& domain_size);
#endif


Real wedge_volume(const Array<dim3,AMREX_SPACEDIM>& elt1,
                  const Array<dim3,AMREX_SPACEDIM>& elt2, 
				  const Vector<int>& is_per_dim, 
				  const Array<Real,AMREX_SPACEDIM>& domain_size);
// -----------------------------------------------------------------------------
// Generic wedge_volume_int for 2D quadrilateral or 3D wedge
//   - DIM = 2 → quadrilateral (A,B in first row; C,D in second row)
//   - DIM = 3 → wedge (A,B,C in first row; D,E,F in second row)
// -----------------------------------------------------------------------------
Real wedge_volume_int(const Array<dim3,AMREX_SPACEDIM>& elt1,
                      const dim3& val1,
                      const Array<dim3,AMREX_SPACEDIM>& elt2,
                      const dim3& val2, 
					  const Vector<int>& is_per_dim, 
					  const Array<Real,AMREX_SPACEDIM>& domain_size);
