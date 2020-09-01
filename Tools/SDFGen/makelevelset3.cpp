#include "makelevelset3.h"
#include <cmath>
#include <AMReX_REAL.H>
#include <utility>
// find distance x0 is from segment x1-x2

//using namespace amrex;

static amrex::Real point_segment_distance(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2)
{
    Vec3r dx(x2-x1);
    amrex::Real m2=mag2(dx);
   // find parameter value of closest point on segment
    amrex::Real s12=(amrex::Real)(dot(x2-x0, dx)/m2);
    if(s12<0){
        s12=0;
    }else if(s12>1){
        s12=1;
    }
   // and find the distance
   return dist(x0, s12*x1+(1-s12)*x2);
}

// find distance x0 is from triangle x1-x2-x3
//static amrex::Real point_triangle_distance(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2, const Vec3r &x3, const Vec3r &normal)
static std::pair<bool,amrex::Real> point_triangle_distance(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2, const Vec3r &x3, const Vec3r &normal,const amrex::Real dx)
{
   // first find barycentric coordinates of closest point on infinite plane
   Vec3r x13(x1-x3), x23(x2-x3), x03(x0-x3);
   amrex::Real m13=mag2(x13), m23=mag2(x23), d=dot(x13,x23);
   amrex::Real invdet=1.f/max(m13*m23-d*d,1e-30);
   amrex::Real a=dot(x13,x03), b=dot(x23,x03);
   // the barycentric coordinates themselves
   amrex::Real w23=invdet*(m23*a-d*b);
   amrex::Real w31=invdet*(m13*b-d*a);
   amrex::Real w12=1-w23-w31;
   
   bool judge=true;
   
   amrex::Real D;
  
 
   Vec3r center;
   //Vec3f localdirection(x0-center);
   for(int i=0;i<3;i++)
   {
       center[i]=(x1[i]+x2[i]+x3[i])/3;
   }
   
   
   amrex::Real dotProduct = ( x0[0]-center[0])*normal[0]+( x0[1]-center[1])*normal[1]+( x0[2]-center[2])*normal[2];
   amrex::Real sign;
  
   //std::cout<<" normal square = "<< normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]<<std::endl;
    
   if(dotProduct>0)
   {
       sign =1.0;
   }
   else
   {
       sign = -1.0;
   }


   if(w23>=0 && w31>=0 && w12>=0)
   { // if we're inside the triangle
      
       D =  dist(x0, w23*x1+w31*x2+w12*x3); 
       return std::make_pair(judge, sign * D); 
        
   }else{ // we have to clamp to one of the edges
      
      if(w23>0)
      { // this rules out edge 2-3 for us
          D =  min(point_segment_distance(x0,x1,x2), point_segment_distance(x0,x1,x3));
        
      }
      else if(w31>0)
      { // this rules out edge 1-3
          D =  min(point_segment_distance(x0,x1,x2), point_segment_distance(x0,x2,x3));
      }
      else
      { // w12 must be >0, ruling out edge 1-2
          D =  min(point_segment_distance(x0,x1,x3), point_segment_distance(x0,x2,x3));
      }
      

      if(D > 0.005 * dx && std::abs(dotProduct) < 0.005 * dx)
      {
           judge = false;
           /*std::cout<<" ========================================"<<std::endl;           
           std::cout<<"parallel points = "<<x0[0]<<"  "<<x0[1]<<"  "<<x0[2]<<"  "<<std::endl;

           std::cout<<"triangle points = "<<x1[0]<<"  "<<x1[1]<<"  "<<x1[2]<<"  "<<std::endl;
           std::cout<<"triangle points = "<<x2[0]<<"  "<<x2[1]<<"  "<<x2[2]<<"  "<<std::endl;
           std::cout<<"triangle points = "<<x3[0]<<"  "<<x3[1]<<"  "<<x3[2]<<"  "<<std::endl;

           
           std::cout<<"dotProduct = "<<dotProduct/dx<<std::endl;
           std::cout<<"dx = "<<dx<<std::endl;
           std::cout<<" ========================================"<<std::endl;
           */


           std::pair<bool,amrex::Real> pd = std::make_pair(judge, sign * D);
           return pd;
      }
      else
      {
           std::pair<bool,amrex::Real> pd = std::make_pair(judge ,sign * D);
           return pd;
      }
   }
   
}

/*//compute normal direction 
Vec3f normal(const Vec3f &x0, const Vec3f &x1, const Vec3f &x2)
{
    Vec3f x12(x1-x0),x23(x2-x1);
    Vec3f normal;

    normal[0] = x12[1]*x23[2]-x12[2]*x23[1]; 
    normal[1] = x12[2]*x23[0]-x12[0]*x23[2];
    normal[2] = x12[0]*x23[1]-x12[1]*x23[0];

    return normal;
}
*/




static void check_neighbour(const std::vector<Vec3ui> &tri, const std::vector<Vec3r> &x,
                            Array3r &phi, Array3i &closest_tri,
                            const Vec3r &gx, int i0, int j0, int k0, int i1, int j1, int k1, const std::vector<Vec3r> &normal,const amrex::Real dx)
{
   if(closest_tri(i1,j1,k1)>=0){
      int s = closest_tri(i1,j1,k1);
      unsigned int p, q, r; assign(tri[s], p, q, r);
      
      std::pair<bool,amrex::Real> d = point_triangle_distance(gx, x[p], x[q], x[r],normal[s],dx);
      
      if(std::abs(d.second)<std::abs(phi(i0,j0,k0)) && d.first == true){
         phi(i0,j0,k0)=d.second;
         closest_tri(i0,j0,k0)=closest_tri(i1,j1,k1);
      }
   }
}

static void sweep(const std::vector<Vec3ui> &tri, const std::vector<Vec3r> &x,
                  Array3r &phi, Array3i &closest_tri, const Vec3r &origin, amrex::Real dx,
                  int di, int dj, int dk,const std::vector<Vec3r> &normal)
{
   int i0, i1;
   if(di>0){ i0=1; i1=phi.ni; }
   else{ i0=phi.ni-2; i1=-1; }
   int j0, j1;
   if(dj>0){ j0=1; j1=phi.nj; }
   else{ j0=phi.nj-2; j1=-1; }
   int k0, k1;
   if(dk>0){ k0=1; k1=phi.nk; }
   else{ k0=phi.nk-2; k1=-1; }
   for(int k=k0; k!=k1; k+=dk) for(int j=j0; j!=j1; j+=dj) for(int i=i0; i!=i1; i+=di){
      Vec3r gx(i*dx+origin[0], j*dx+origin[1], k*dx+origin[2]);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i-di, j,    k,normal,dx);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i,    j-dj, k,normal,dx);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i-di, j-dj, k,normal,dx);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i,    j,    k-dk,normal,dx);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i-di, j,    k-dk,normal,dx);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i,    j-dj, k-dk,normal,dx);
      check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i-di, j-dj, k-dk,normal,dx);
   }
}

// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
static int orientation(amrex::Real x1, amrex::Real y1, amrex::Real x2, amrex::Real y2, amrex::Real &twice_signed_area)
{
   twice_signed_area=y1*x2-x1*y2;
   if(twice_signed_area>0) return 1;
   else if(twice_signed_area<0) return -1;
   else if(y2>y1) return 1;
   else if(y2<y1) return -1;
   else if(x1>x2) return 1;
   else if(x1<x2) return -1;
   else return 0; // only true when x1==x2 and y1==y2
}

// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
// if true is returned, the barycentric coordinates are set in a,b,c.
static bool point_in_triangle_2d(amrex::Real x0, amrex::Real y0, 
                                 amrex::Real x1, amrex::Real y1, amrex::Real x2, amrex::Real y2, amrex::Real x3, amrex::Real y3,
                                 amrex::Real& a, amrex::Real& b, amrex::Real& c)
{
   x1-=x0; x2-=x0; x3-=x0;
   y1-=y0; y2-=y0; y3-=y0;
   int signa=orientation(x2, y2, x3, y3, a);
   if(signa==0) return false;
   int signb=orientation(x3, y3, x1, y1, b);
   if(signb!=signa) return false;
   int signc=orientation(x1, y1, x2, y2, c);
   if(signc!=signa) return false;
   amrex::Real sum=a+b+c;
   assert(sum!=0); // if the SOS signs match and are nonkero, there's no way all of a, b, and c are zero.
   a/=sum;
   b/=sum;
   c/=sum;
   return true;
}

void make_level_set3(const std::vector<Vec3ui> &tri, const std::vector<Vec3r> &x, const std::vector<Vec3r> &normal,
                     const Vec3r &origin, amrex::Real dx, int ni, int nj, int nk,
                     Array3r &phi, amrex::Real dmax, const int exact_band)
{
   phi.resize(ni, nj, nk);
   phi.assign(dmax); // upper bound on distance

   Array3i closest_tri(ni, nj, nk, -1);
   Array3i intersection_count(ni, nj, nk, 0); // intersection_count(i,j,k) is # of tri intersections in (i-1,i]x{j}x{k}
   // we begin by initializing distances near the mesh, and figuring out intersection counts
   Vec3r ijkmin, ijkmax;
   
   for(unsigned int t=0; t<tri.size(); ++t){
     unsigned int p, q, r; assign(tri[t], p, q, r);
     // coordinates in grid to high precision
      amrex::Real fip=((amrex::Real)x[p][0]-origin[0])/dx, fjp=((amrex::Real)x[p][1]-origin[1])/dx, fkp=((amrex::Real)x[p][2]-origin[2])/dx;
      amrex::Real fiq=((amrex::Real)x[q][0]-origin[0])/dx, fjq=((amrex::Real)x[q][1]-origin[1])/dx, fkq=((amrex::Real)x[q][2]-origin[2])/dx;
      amrex::Real fir=((amrex::Real)x[r][0]-origin[0])/dx, fjr=((amrex::Real)x[r][1]-origin[1])/dx, fkr=((amrex::Real)x[r][2]-origin[2])/dx;
      // do distances nearby
      int i0=clamp(int(min(fip,fiq,fir))-exact_band, 0, ni-1), i1=clamp(int(max(fip,fiq,fir))+exact_band+1, 0, ni-1);
      int j0=clamp(int(min(fjp,fjq,fjr))-exact_band, 0, nj-1), j1=clamp(int(max(fjp,fjq,fjr))+exact_band+1, 0, nj-1);
      int k0=clamp(int(min(fkp,fkq,fkr))-exact_band, 0, nk-1), k1=clamp(int(max(fkp,fkq,fkr))+exact_band+1, 0, nk-1);
      for(int k=k0; k<=k1; ++k) for(int j=j0; j<=j1; ++j) for(int i=i0; i<=i1; ++i){
         Vec3r gx(i*dx+origin[0], j*dx+origin[1], k*dx+origin[2]);
         std::pair<bool,amrex::Real> d=point_triangle_distance(gx, x[p], x[q], x[r], normal[t],dx);
	 //std::cout<<"d="<<d<<std::endl;
         //std::cout<<"dx="<<dx<<std::endl;
         if(std::abs(d.second)<std::abs(phi(i,j,k)) && d.first == true){
   //              if(d==0)
   //              {
   //                   std::cout<<"d="<<d<<std::endl;
   //		 }
                 phi(i,j,k)=d.second;
		 closest_tri(i,j,k)=t;
	 }
      }
      // and do intersection counts
      j0=clamp((int)std::ceil(min(fjp,fjq,fjr)), 0, nj-1);
      j1=clamp((int)std::floor(max(fjp,fjq,fjr)), 0, nj-1);
      k0=clamp((int)std::ceil(min(fkp,fkq,fkr)), 0, nk-1);
      k1=clamp((int)std::floor(max(fkp,fkq,fkr)), 0, nk-1);
      for(int k=k0; k<=k1; ++k) for(int j=j0; j<=j1; ++j){
	      amrex::Real a, b, c;
	      if(point_in_triangle_2d(j, k, fjp, fkp, fjq, fkq, fjr, fkr, a, b, c)){
		      amrex::Real fi=a*fip+b*fiq+c*fir; // intersection i coordinate
		      int i_interval = int(std::ceil(fi)); // intersection is in (i_interval-1,i_interval]
		      if(i_interval<0) ++intersection_count(0, j, k); // we enlarge the first interval to include everything to the -x direction
		      else if(i_interval<ni) ++intersection_count(i_interval,j,k);
		      // we ignore intersections that are beyond the +x side of the grid
	      }
      }
   }
   // and now we fill in the rest of the distances with fast sweeping
   for(unsigned int pass=0; pass<2; ++pass){
      sweep(tri, x, phi, closest_tri, origin, dx, +1, +1, +1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, -1, -1, -1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, +1, +1, -1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, -1, -1, +1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, +1, -1, +1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, -1, +1, -1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, +1, -1, -1,normal);
      sweep(tri, x, phi, closest_tri, origin, dx, -1, +1, +1,normal);
   }
#if 0
   // then figure out signs (inside/outside) from intersection counts
   for(int k=0; k<nk; ++k) for(int j=0; j<nj; ++j){
      int total_count=0;
      for(int i=0; i<ni; ++i){
         total_count+=intersection_count(i,j,k);
         if(total_count%2==1){ // if parity of intersections so far is odd,
            phi(i,j,k)=-phi(i,j,k); // we are inside the mesh
         }
      }
   }
#endif
}

