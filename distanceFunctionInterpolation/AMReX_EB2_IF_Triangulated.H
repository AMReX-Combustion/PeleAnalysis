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

namespace amrex { namespace EB2 {

struct pointToElement
    {
    public: 
        Real physLocation[BL_SPACEDIM];
        long elementIndex;
        int pointIndex;    
        
        int mergeTag;

        pointToElement(std::vector<Real> pt,long ei, int pi):
        elementIndex(ei),pointIndex(pi),mergeTag(0)
        {
            physLocation[0]=pt[0];
            physLocation[1]=pt[1];
#if BL_SPACEDIM==3            
            physLocation[2]=pt[2];
#endif
        }



    bool operator <(const pointToElement& rhs)const
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
        else
        {
            return this->physLocation[0]<rhs.physLocation[0];
        }       

    }



    bool operator >(const pointToElement& rhs)const 
    {
        if(this->physLocation[0]==rhs.physLocation[0])
        {

#if BL_SPACEDIM==3            
            
            if(this->physLocation[1]==rhs.physLocation[1])
            {
                return this->physLocation[2]>rhs.physLocation[2];
            }
            
            else
#endif
            {
                return this->physLocation[1]>rhs.physLocation[1];
            }
        }
        else
        {
            return this->physLocation[0]>rhs.physLocation[0];
        }    
    }

    bool operator ==(const pointToElement& rhs)const
    {
        return (this->physLocation[0]==rhs.physLocation[0] && this->physLocation[1]==rhs.physLocation[1] 

#if BL_SPACEDIM==3
           
               && this->physLocation[2]==rhs.physLocation[2]

#endif
               );

    }

    bool operator !=(const pointToElement& rhs)const
    {
        return (this->physLocation[0]!=rhs.physLocation[0] || this->physLocation[1]!=rhs.physLocation[1] 

#if BL_SPACEDIM==3

                    || this->physLocation[2]!=rhs.physLocation[2]

#endif

               );
    }

 

};

class TriangulatedIF
//      : public GPUable
//      : public IF_base....etc...

{
public:
    //Triangulated Implicit Function: should be built with triangulated surface file
    //extendable for liquid filled surface
    
    //----------Constructors--------------------------
    
       //----Construct from file----------- 
   // TriangulatedIF(const char* nameOfFile, const char* type);

    TriangulatedIF(const char* isoFile);

    //---------Destructors----------------------------
    ~TriangulatedIF(){};

//protected:
    
    
   //---------Member Functions-----------------------


//Data loading 

       //load data from non-watertight format
    void loadData (const char* nameOfFile , std::vector<std::vector<Real> >& temp_surface);              
       
       //re-organize and make it water-tight
    void reOrganize
    (
          const std::vector<std::vector<Real> >& temp_surface
    );

    void mergeVertex ( std::vector<pointToElement>& temp_ptlist );
    
    Real distance
    (
        const Real (&p1)[BL_SPACEDIM],
        const Real (&p2)[BL_SPACEDIM] 
    )const;

//Distance Function building

    void buildDistance(char *isoFild);
    
    
    //---------Variables-----------------------------
    
    //--------normal 
    std::vector<std::vector<Real> > normalList;
  
    //--------point index
    std::vector<std::vector<long> > faceList; 

   //--------Point Locations
    std::vector<std::vector<Real> > vertList; 

    double eps;

    MultiFab Dfab;
 
    Geometry geom;
    
    amrex::distanceFunctionInterpolation distanceInterpolation; 
};
   

}}


#endif
