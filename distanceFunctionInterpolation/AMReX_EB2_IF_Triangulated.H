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
      TriangulatedIF(const std::string& isoFile);

      ~TriangulatedIF();

      //load data from non-watertight format
      void loadData (const std::string&               nameOfFile,
                     std::vector<std::vector<Real> >& temp_surface);              
       
      //re-organize and make it water-tight
      void reOrganize (const std::vector<std::vector<Real> >& temp_surface);

      void mergeVertex (std::vector<pointToElement>& temp_ptlist);
    
      Real distance (const Array<Real,AMREX_SPACEDIM> &p1,
                     const Array<Real,AMREX_SPACEDIM> &p2 )const;

      void buildDistance(const std::string& IsoFile);
    
      std::vector<std::vector<Real> > normalList;
      std::vector<std::vector<int> > faceList; 
      std::vector<std::vector<Real> > vertList; 
      double eps;
      MultiFab distanceMF; 
      Geometry geom;
      distanceFunctionInterpolation* distanceInterpolation_; 
    };
  }
}


#endif
