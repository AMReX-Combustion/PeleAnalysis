#ifndef AMREX_EB2_IF_TRIANGULATED_H_
#define AMREX_EB2_IF_TRIANGULATED_H_

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Base.H>
#include <vector>
#include <string.h>
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
      TriangulatedIF(const std::string& isoFile, const std::string& fileType);
        
      AMREX_GPU_HOST_DEVICE inline
      Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept {
        return Distance(AMREX_D_DECL(x,y,z));
      }

      inline Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
      }  

      void finalize(const Geometry& geom, const BoxArray& grids, const DistributionMapping& dm);     
 
      const DistributionMapping& dm()
      {
          return dm_;
      }
    
      const Geometry& geom()
      {
          return geom_; 
      }

      const BoxArray& grids()
      {
          return grids_;
      }

    protected:
      //load data from non-watertight format
      void loadData(const std::string& isoFile, const std::string& fileType) noexcept;              
      void loadData_mef(const std::string& isoFile) noexcept;      
      void loadData_stl_ascii (const std::string& isoFile) noexcept;
      void loadData_stl_binary(const std::string& isoFile) noexcept;

      void makeWatertight (const std::vector<std::vector<Real> >& temp_surface) noexcept;

      void mergeVertex (std::vector<pointToElement>& temp_ptlist) noexcept;
    
      void buildDistance() noexcept;

      Real Distance(AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
        {
          const Real* dxInv = geom_.InvCellSize(); 
          IntVect iv(AMREX_D_DECL(floor(x * dxInv[0]),
                                  floor(y * dxInv[1]),
                                  floor(z * dxInv[2])));

          Real distance;
          for (MFIter mfi(distanceMF_); mfi.isValid(); ++mfi)
          {
            if (mfi.validbox().contains(iv))
            {            
              const auto& d=distanceMF_.array(mfi);
              const Box& vbox = mfi.tilebox();

              const Real n[AMREX_SPACEDIM] = {AMREX_D_DECL(x*dxInv[0] - iv[0],
                                                           y*dxInv[1] - iv[1],
                                                           z*dxInv[2] - iv[2])};

#if AMREX_SPACEDIM==2
              distance =/*0 0*/ d(iv[0]  ,iv[1]  ,iv[2]  ) * (1 - n[0]) * (1 - n[1]) 
                       +/*1 0*/ d(iv[0]+1,iv[1]  ,iv[2]  ) * n[0]       * (1 - n[1]) 
                       +/*0 1*/ d(iv[0]  ,iv[1]+1,iv[2]  ) * (1 - n[0]) * n[1]
                       +/*1 1*/ d(iv[0]  ,iv[1]  ,iv[2]+1) * n[0]       * n[1];
#else
              distance =/*0 0 0*/ d(iv[0]  ,iv[1]  ,iv[2]  ) * (1 - n[0]) * (1 - n[1]) * (1 - n[2])
                       +/*1 0 0*/ d(iv[0]+1,iv[1]  ,iv[2]  ) * n[0]       * (1 - n[1]) * (1 - n[2])
                       +/*0 1 0*/ d(iv[0]  ,iv[1]+1,iv[2]  ) * (1 - n[0]) * n[1]       * (1 - n[2])
                       +/*0 0 1*/ d(iv[0]  ,iv[1]  ,iv[2]+1) * (1 - n[0]) * (1 - n[1]) * n[2]
                       +/*1 1 0*/ d(iv[0]+1,iv[1]+1,iv[2]  ) * n[0]       *     n[1]   * (1 - n[2])
                       +/*1 0 1*/ d(iv[0]+1,iv[1]  ,iv[2]+1) * n[0]       * (1 - n[1]) * n[2]
                       +/*0 1 1*/ d(iv[0]  ,iv[1]+1,iv[2]+1) * (1 - n[0]) * n[1]       * n[2]
                       +/*1 1 1*/ d(iv[0]+1,iv[1]+1,iv[2]+1) * n[0]       * n[1]       * n[2];
#endif
              break;               
            }        
          }
          return distance;
        }
      
      
      //-------------members--------------------------
      std::vector<Vec3r> normalList;
   
      std::vector<Vec3ui> faceList; 
   
      std::vector<Vec3r> vertList; 
   
      double eps;

      MultiFab distanceMF_; 

      DistributionMapping dm_; 

      Geometry geom_;

      BoxArray grids_;
    };
  }
}


#endif
