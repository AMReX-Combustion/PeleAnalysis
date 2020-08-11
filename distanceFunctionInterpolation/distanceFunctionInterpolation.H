class distanceFunctionInterpolation
{
public:

//--------------------------------------constructor-----------------------------------//
distanceFunctionInterpolation(MultiFab& dfab, Geometry& geom)
: distanceFab_(dfab),geom_(geom){}


//----------------------------------- member function-------------------------------- //

Real distance(XDim3 location)
{
   // boxArray ba=distanceFab.bxs;
    Real distance;
 //   IntVect small,big;

    IntVect iv;

    dx = geom.CellSize(); 

    for(i=0;i<AMREX_DIM;i++)
    {
        iv[i]=floor(location[i]/dx[i]);
    }

    for(MFIter mfi(distanceFab_.bxs,distanceFab_.dm, true);mfi.isValid();++mfi)
    {
        if(mfi.contains(location))
        {
            IntVect small(mfi.smallend),big(mfi.bigend);
            Real ksai,yita,zita,side;
            
            side = (Real)big[2] - (Real)small[2];
            ksai = ( location[0]-(Real)small[0] )/side;
            yita = ( location[1]-(Real)small[1] )/side;
            zita = ( location[2]-(Real)small[2] )/side;

            distance =/*0 0 0*/ a1 * (1 - ksai) * (1 - yita) * (1 - zita)
                     +/*1 0 0*/ a2 * ksai       * (1 - yita) * (1 - zita)
                     +/*0 1 0*/ a3 * (1 - ksai) * yita       * (1 - zita)
                     +/*0 0 1*/ a4 * (1 - ksai) * (1 - yita) * zita
                     +/*1 1 0*/ a5 * ksai       *     yita   * (1 - zita)
                     +/*1 0 1*/ a6 * ksai       * (1 - yita) * zita
                     +/*0 1 1*/ a7 * (1 - ksai) * yita       * zita
                     +/*1 1 1*/ a8 * ksai       * yita       * zita;
           

            break;               
                      
        }
    }
    return distance;
}




//----------------------------------- protected member function --------------------- //
MultiFab& distanceMap()
{
    return distanceFab_;
}
Geometry& geom()
{
    return geom_;
}

protected:

MultiFab& distanceFab_;

Geometry& geom_;
}
