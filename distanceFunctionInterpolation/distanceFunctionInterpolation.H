#ifndef AMREX_EB2_IF_DISTANCEINTERPOLATION_H_
#define AMREX_EB2_IF_DISTANCEINTERPOLATION_H_

#include <AMReX_EB2_IF_Triangulated.H>
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>


#include <AMReX_ParmParse.H>

#include <AMReX_MultiFab.H>

namespace amrex{ namespace EB2 {
class distanceFunctionInterpolation
{
public:

//--------------------------------------constructor-----------------------------------//
distanceFunctionInterpolation(MultiFab& dfab, const Geometry& geom)
: distanceFab_(dfab),geom_(geom){}

//----------------------------------- member function-------------------------------- //



Real distance(Real* location)
{
   // boxArray ba=distanceFab.bxs;
    Real distance;
 //   IntVect small,big;

    IntVect iv;

    const Real* dx = geom_.CellSize(); 

    for(int i=0;i<BL_SPACEDIM;i++)
    {
        iv[i]=floor(location[i]/dx[i]);
    }

    

    for(MFIter mfi(distanceFab_.boxarray,distanceFab_.distributionMap);mfi.isValid();++mfi)
    {
        if(mfi.validbox().contains(iv))
        {
            
            const auto& d=distanceFab_.array(mfi);

            const Box& vbox = distanceFab_.boxarray[mfi.index()];

            const int* lowIndex  = vbox.loVect();

            const int* highIndex = vbox.hiVect(); 

            Real ksai,yita,zita;
            
            ksai = location[0]/dx[0]-iv[0];
            yita = location[1]/dx[1]-iv[1];
            zita = location[2]/dx[2]-iv[2];

            distance =/*0 0 0*/ d(iv[0]  ,iv[1]  ,iv[2]  ) * (1 - ksai) * (1 - yita) * (1 - zita)
                     +/*1 0 0*/ d(iv[0]+1,iv[1]  ,iv[2]  ) * ksai       * (1 - yita) * (1 - zita)
                     +/*0 1 0*/ d(iv[0]  ,iv[1]+1,iv[2]  ) * (1 - ksai) * yita       * (1 - zita)
                     +/*0 0 1*/ d(iv[0]  ,iv[1]  ,iv[2]+1) * (1 - ksai) * (1 - yita) * zita
                     +/*1 1 0*/ d(iv[0]+1,iv[1]+1,iv[2]  ) * ksai       *     yita   * (1 - zita)
                     +/*1 0 1*/ d(iv[0]+1,iv[1]  ,iv[2]+1) * ksai       * (1 - yita) * zita
                     +/*0 1 1*/ d(iv[0]  ,iv[1]+1,iv[2]+1) * (1 - ksai) * yita       * zita
                     +/*1 1 1*/ d(iv[0]+1,iv[1]+1,iv[2]+1) * ksai       * yita       * zita;
           

            break;               
                      
        }
    }
    return distance;
}




//----------------------------------- protected member function --------------------- //
const MultiFab& distanceMap()
{
    return distanceFab_;
}
const Geometry& geom()
{
    return geom_;
}

protected:

const MultiFab& distanceFab_;

const Geometry& geom_;
};

}}
#endif
