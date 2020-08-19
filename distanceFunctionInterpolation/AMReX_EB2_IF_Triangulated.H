#ifndef AMREX_EB2_IF_TRIANGULATED_H_
#define AMREX_EB2_IF_TRIANGULATED_H_

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Base.H>
#include <vector>
#include <string.h>
#include "distanceFunctionInterpolation.H"
#include <AMReX_MultiFab.H>
#include <vec.h>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex
{
  namespace EB2
  {
    struct pointToElement
    {
    public: 
      Array<Real,AMREX_SPACEDIM> physLocation;
      int elementIndex;
      int pointIndex;
      int mergeTag;

      pointToElement(std::vector<Real> pt,int ei, int pi)
        : elementIndex(ei),
          pointIndex(pi),
          mergeTag(0)
        {
          physLocation = {D_DECL(pt[0],pt[1],pt[2])};
        }

      bool operator <(const pointToElement& rhs) const
        {
          if(this->physLocation[0]==rhs.physLocation[0])
          {
#if BL_SPACEDIM==3                        
            if(this->physLocation[1]==rhs.physLocation[1])
            {
              return this->physLocation[2]<rhs.physLocation[2];
            }
            else
#endif
            {
              return this->physLocation[1]<rhs.physLocation[1];
            }
          }
          return this->physLocation[0]<rhs.physLocation[0];
        }
    };

    class TriangulatedIF
    {
    public:
      //----Construct from file----------- 
      TriangulatedIF(/*const Geometry& geom, const BoxArray& grids,*/ const std::string& isoFile, const std::string& fileType);

      ~TriangulatedIF();

      //load data from non-watertight format
      void loadData(const std::string& isoFile, const std::string& fileType);

      void loadData_stl_ascii (const std::string& isoFile);
              
      void loadData_mef(const std::string& isoFile); 
     
      void loadData_stl_binary(const std::string& isoFile);

      //re-organize and make it water-tight
      void reOrganize (const std::vector<std::vector<Real> >& temp_surface);

      void mergeVertex (std::vector<pointToElement>& temp_ptlist);
    
      Real distance (const Array<Real,AMREX_SPACEDIM> &p1,
                     const Array<Real,AMREX_SPACEDIM> &p2 )const;

      void buildDistance(/*const std::string& isoFile*/);
    
      void finalize(const Geometry& geom, const BoxArray& grids, const DistributionMapping& dm);     
 
      const distanceFunctionInterpolation& distanceInterpolation()
      {
          return *distanceInterpolation_;
      }
      
      const DistributionMapping& dm()
      {
          return *dm_;
      }
    
      const Geometry& geom()
      {
          return *geom_; 
      }

      const BoxArray& grids()
      {
          return *grids_;
      }
//-------------members--------------------------
       std::vector<Vec3r> normalList;
   
      std::vector<Vec3ui> faceList; 
   
      std::vector<Vec3r> vertList; 
   
      double eps;

      MultiFab distanceMF; 

      const DistributionMapping* dm_; 

      const Geometry* geom_;

      const BoxArray* grids_;

      distanceFunctionInterpolation* distanceInterpolation_;


      
    };
  }
}


#endif
