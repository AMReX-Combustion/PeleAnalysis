//#ifndef AMREX_EB2_IF_SPHERE_H_
//#define AMREX_EB2_IF_SPHERE_H_
#include <AMReX_EB2_IF_Triangulated.H>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <cmath>
#include <algorithm>
// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

//this class is designed to call std::sort to sort the generated pointlist
//once sorted, merge vertices 



//----------Constructors--------------------------
    
/*       //----Construct from file----------- 
    TriangulatedIF::TriangulatedIF(const char* nameOfFile, const char* type)
    {
         std::vector<std::vector<Real> > temp_surface;

         TriangulatedIF::loadData(nameOfFile,temp_surface);
          
         TriangulatedIF::reOrganize(temp_surface); 
     
     //  call makelevelset3.cpp and return is as a multifab
         
         
    }
*/
    TriangulatedIF::TriangulatedIF(const char* isoFile)
    {
    //     std::vector<std::vector<Real> > normalList;

    //     std::vector<std::vector<long> > faceList;

    //     std::vector<std::vector<Real> > vertList; 

     //    Geometry geom;
  
     //    MultiFab Dfab;

         TriangulatedIF::buildDistance(isoFile);
         
         distanceInterpolation_=new distanceFunctionInterpolation(Dfab,geom);      
    }

    TriangulatedIF::~TriangulatedIF()
    {
        delete distanceInterpolation_;
    }
    //---------Protected Member Functions-----------------------
    void TriangulatedIF::loadData
    (
          const char* nameOfFile,
          std::vector<std::vector<Real> >& temp_surface 
          
    ) 
    {
          //designed for stl file, if needed, use factory for more extensions
          
          //load data
          FILE *fp=NULL;
          
          fp=fopen(nameOfFile,"r");
         
          double num1,num2,num3;
          
          int max_line=512;

          char dump[max_line];
          //Title line, discarded
          fgets(dump,max_line,(FILE*)fp);
          
          for(int i=0;;++i)
          {
               
               fscanf(fp,"%s",dump);
               
               if(strcmp(dump,"facet")!=0)
               {
                   break;    
               }
               fscanf(fp,"%s",dump); 
              
               this->normalList.push_back(std::vector<Real>());
 
               fscanf(fp,"%lf %lf %lf\n",&num1,&num2,&num3);
                
               this->normalList[i].push_back( (Real)num1 );
               this->normalList[i].push_back( (Real)num2 );
               this->normalList[i].push_back( (Real)num3 );         
              
               fgets(dump,max_line,(FILE*)fp);
               
               for(int j=0;j<3;++j)
               {
                   fscanf(fp,"%s %lf %lf %lf\n",dump,&num1,&num2,&num3);
                  
                   temp_surface.push_back(std::vector<Real>());
           

                   temp_surface[3*i+j].push_back( (Real)num1 );
                   temp_surface[3*i+j].push_back( (Real)num2 );
                   temp_surface[3*i+j].push_back( (Real)num3 );
               }

               fgets(dump,max_line,(FILE*)fp);
               fgets(dump,max_line,(FILE*)fp);
          }
          fclose(fp); 
          
    }
                
// do simplification and make data water-tight             
    void TriangulatedIF::reOrganize
    (
          const std::vector<std::vector<Real> >& temp_surface
    )
    {
  
        std::vector< pointToElement > temp_ptlist;

        for(int i=0;i<this->normalList.size();++i)
        {
            faceList.push_back(std::vector<int>() );           
            for(int j=0;j<3;++j)
            {
                faceList[i].push_back(3*i+j);

                temp_ptlist.push_back(pointToElement(temp_surface[3*i+j],i,j) );                        
            }
        }
        std::sort(temp_ptlist.begin(),temp_ptlist.end());
     
        mergeVertex(temp_ptlist);
    }
    
    // use a stupid algorithm here to merge vertices, acceleration may applied
    
    void TriangulatedIF::mergeVertex(std::vector<pointToElement>& temp_ptlist)
    {
        int i=0,j=0,indexCounter=0,stepCounter;
        
        int mergeNumber=0;
        
        for(i=0;i<temp_ptlist.size();++i)
        {
            stepCounter=0;

            mergeNumber=0;
      
            vertList.push_back(std::vector<Real>());
            
            vertList[indexCounter].push_back(temp_ptlist[i].physLocation[0]);
            vertList[indexCounter].push_back(temp_ptlist[i].physLocation[1]);
#if BL_SPACEDIM==3    
          
            vertList[indexCounter].push_back(temp_ptlist[i].physLocation[2]);

#endif            
            faceList[temp_ptlist[i].elementIndex][temp_ptlist[i].pointIndex] = indexCounter;
           
            for(j=i+1;j<temp_ptlist.size();++j)
            {
                //try to reduce the amount of point should go through, if xj>xi+eps and yj>yi+eps, distance(i,j)>eps
       /*         if(temp_ptlist[j].physLocation[0]>temp_ptlist[i].physLocation[0]+eps && temp_ptlist[j].physLocation[1]>temp_ptlist[i].physLocation[1]+eps)
                {
                    if(j==i+1)
                    {
                        std::cout<<"give stepCounter a value"<<std::endl;
                        stepCounter=1;
                        indexCounter++;
                    }
                    break;
                }*/
                if(distance(temp_ptlist[i].physLocation,temp_ptlist[j].physLocation)<eps)
                {
                    //merge vertex
                    faceList[temp_ptlist[j].elementIndex][temp_ptlist[j].pointIndex] = indexCounter;
                    
                    temp_ptlist[j].mergeTag==1;
                    
                }
                else
                { 
                    if(temp_ptlist[j].mergeTag==0)
                    {
                        if(mergeNumber==0)
                        {
                            stepCounter=j-i;
                            indexCounter++;
                            mergeNumber++;
                        }
                    }
                    else
                    {}
                }
            i=stepCounter-1;
            }
            
        }
          
    }    
    
    
    //return distance between point p1 and point p2
    Real TriangulatedIF::distance
    (
        const Real (&p1)[BL_SPACEDIM],
        const Real (&p2)[BL_SPACEDIM]     
    )const
    {
        return std::sqrt(pow((p1[0]-p2[0]),2)+pow((p1[1]-p2[1]),2)
        
#if BL_SPACEDIM==3     
       
                +pow((p1[2]-p2[2]),2)

#endif

                   );
    }
    
/*    std::vector<std::vector<Real> >& TriangulatedIF::normal()
    {
        return &normal;
    } 

    std::vector<std::vector<long> >& TriangulatedIF::surface()
    {
        return &surface;
    }

    std::vector<std::vector<Real> >& TriangulatedIF::pointList()
    {
        return &pointList;
    }
*/ 
}}//endnamespace EB2 & amrex
