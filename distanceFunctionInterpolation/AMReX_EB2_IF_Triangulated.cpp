//#ifndef AMREX_EB2_IF_SPHERE_H_
//#define AMREX_EB2_IF_SPHERE_H_
#include <AMReX_EB2_IF_Triangulated.H>
#include <stdio.h>
#include <vector>
#include <string>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include "makelevelset3.h"

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex
{
    //file type table,extendable 
    enum FILETYPE
    {
        unknown    = 0,
        mef        = 1,
        stl_ascii  = 2,
        stl_binary = 3
    };
    static std::map<std::string,FILETYPE> typeMap = {{"mef", mef}, {"stl_ascii", stl_ascii}, {"stl_binary", stl_binary}};

    Vec3r normal(const Vec3r &x0, const Vec3r &x1, const Vec3r &x2)
    {
        Vec3r x12(x1-x0),x23(x2-x1);
        Vec3r normal;

        normal[0] = x12[1]*x23[2]-x12[2]*x23[1];
        normal[1] = x12[2]*x23[0]-x12[0]*x23[2];
        normal[2] = x12[0]*x23[1]-x12[1]*x23[0];

        float r = std::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
        normal[0] /= r;
        normal[1] /= r;
        normal[2] /= r;
        return normal;
    }

    static amrex::Real distance (const Array<amrex::Real,AMREX_SPACEDIM> &p1,
                                 const Array<amrex::Real,AMREX_SPACEDIM> &p2)
    {
        return std::sqrt(AMREX_D_TERM(pow((p1[0]-p2[0]),2),+pow((p1[1]-p2[1]),2),+pow((p1[2]-p2[2]),2)));
    }

    static
    std::vector<std::string> parseVarNames(std::istream& is)
    {
        std::string line;
        std::getline(is,line);
        return amrex::Tokenize(line,std::string(", "));
    }

    static std::string parseTitle(std::istream& is)
    {
        std::string line;
        std::getline(is,line);
        return line;
    }

    namespace EB2
    {
        Real TriangulatedIF::eps = 1e-50; // What should this be?
        
        TriangulatedIF::TriangulatedIF(const std::string& isoFile,
                                       const std::string& fileType)
        {
            TriangulatedIF::loadData(isoFile,fileType); 
        }

        //---------Protected Member Functions-----------------------
        void
        TriangulatedIF::finalize(const Geometry&            geom,
                                 const BoxArray&            grids,
                                 const DistributionMapping& dm) noexcept
        {
            geom_  = geom;  
            grids_ = grids;
            dm_    = dm; 

            int nGrow = 1; // nGrow must equal narrowband half-width
            distanceMF_.define(this->grids(), this->dm(), 1, 1);
           
            Print()<<"start build distance"<<std::endl;

            buildDistance();
        }
    
        void
        TriangulatedIF::loadData (const std::string& isoFile,
                                  const std::string& fileType) noexcept
        {
            std::cout<<"fileType="<<typeMap[fileType]<<std::endl;
            switch(typeMap[fileType])
            {
            case mef:
 
                loadData_mef(isoFile);
                break;
 
            case stl_ascii:
                loadData_stl_ascii(isoFile);
                break;

            case stl_binary:       
                loadData_stl_binary(isoFile);
                break;
            
            default:
                Abort("Unknown surface file format");
            }
        }
        
        void
        TriangulatedIF::loadData_mef (const std::string& isoFile) noexcept
        {
            if (ParallelDescriptor::IOProcessor()) 
            {
                std::cerr << "Reading isoFile... " << isoFile << std::endl;
            }
    
            FArrayBox nodes;
            Vector<int> faceData;
            int nElts;
            int nodesPerElt;
            int nNodes, nCompNodes;
            std::vector<std::string> surfNames;
  
            std::ifstream ifs;
            ifs.open(isoFile.c_str(),std::ios::in|std::ios::binary);

            const std::string title = parseTitle(ifs);
            surfNames = parseVarNames(ifs);
            nCompNodes = surfNames.size();

            ifs >> nElts;
            ifs >> nodesPerElt;

            FArrayBox tnodes;
            tnodes.readFrom(ifs);
            nNodes = tnodes.box().numPts();

            // transpose the data so that the components are 'in the right spot for fab data'
            // NOTE: This is unnecessary, but eases conversion below
            nodes.resize(tnodes.box(),nCompNodes);
            Real** np = new Real*[nCompNodes];
            for (int j=0; j<nCompNodes; ++j)
            {
                np[j] = nodes.dataPtr(j);
            }
    
            Real* ndat = tnodes.dataPtr();
    
            for (int i=0; i<nNodes; ++i)
            {
                for (int j=0; j<nCompNodes; ++j)
                {
                    np[j][i] = ndat[j];
                }
                ndat += nCompNodes;
            }
            delete [] np;
            tnodes.clear();

            faceData.resize(nElts*nodesPerElt,0);
            ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
            ifs.close();

            // Convert surface data into form required for distance function computation
            for (int node=0; node<nNodes; ++node) 
            {
                const IntVect iv(D_DECL(node,0,0));
                vertList.push_back(Vec3r(D_DECL(nodes(iv,0),nodes(iv,1),nodes(iv,2))));
            }
            for (int elt=0; elt<nElts; ++elt) 
            {
                int offset = elt * nodesPerElt;
                faceList.push_back(Vec3ui(D_DECL(faceData[offset]-1,faceData[offset+1]-1,faceData[offset+2]-1)));
                normalList.push_back(normal(vertList[faceData[offset]-1],vertList[faceData[offset+1]-1],vertList[faceData[offset+2]-1]));
            }

            Print() << "Isosurface [format: mef] contains " << nElts << " elements and " << nNodes << " nodes" << std::endl;
        }
        
        void
        TriangulatedIF::loadData_stl_ascii (const std::string& isoFile) noexcept
        {
            std::vector<std::vector<Real> > temp_surface;

            FILE *fp=NULL;
            
            fp = fopen(isoFile.c_str(),"r");
            
            double num1,num2,num3;

            int max_line=512;

            char dump[max_line];
            //Title line, discarded
            //
            fgets(dump,max_line,(FILE*)fp);

            for(int i=0;;++i)
            {
                fscanf(fp,"%s",dump);

    //            std::cout<<"dump="<<dump<<std::endl;

                if(strcmp(dump,"facet")!=0)
                {
                    break;    
                }
                fscanf(fp,"%s",dump); 

                fscanf(fp,"%lf %lf %lf\n",&num1,&num2,&num3);

  //              std::cout<<"normal="<<num1<<"    "<<num2<<"    "<<num3<<std::endl;

                this->normalList.push_back(Vec3r((Real)num1,(Real)num2,(Real)num3) );

                fgets(dump,max_line,(FILE*)fp);

                for(int j=0;j<3;++j)
                {
                    fscanf(fp,"%s %lf %lf %lf\n",dump,&num1,&num2,&num3);

                    temp_surface.push_back(std::vector<Real>());

                    temp_surface[3*i+j].push_back( (Real)num1 );
                    temp_surface[3*i+j].push_back( (Real)num2 );
                    temp_surface[3*i+j].push_back( (Real)num3 );
//                std::cout<<"p"<<j+1<<"="<<num1<<"    "<<num2<<"    "<<num3<<std::endl;
                }

                fgets(dump,max_line,(FILE*)fp);
                fgets(dump,max_line,(FILE*)fp);
            }
            fclose(fp);

            makeWatertight(temp_surface);
        }
   
        void
        TriangulatedIF::loadData_stl_binary (const std::string& isoFile) noexcept
        {
            std::vector<std::vector<Real> > temp_surface;

            FILE *fp = fopen(isoFile.c_str(),"rb");
          
            char dump[80];

            fread(dump,80,1,fp);          
 
            uint32_t faceCount;
       
            fread(&faceCount, sizeof(uint32_t), 1, fp);
               
            float num1,num2,num3;

            char c[2];
            for (unsigned int i = 0; i < faceCount; i++)
            {
                fread(&num1,sizeof(float),1,fp);
                fread(&num2,sizeof(float),1,fp);
                fread(&num3,sizeof(float),1,fp);

                this->normalList.push_back(Vec3r((Real)num1,(Real)num2,(Real)num3) ); 

                for(int j=0;j<3;++j)
                {                   
                    fread(&num1,sizeof(float),1,fp);
                    fread(&num2,sizeof(float),1,fp);
                    fread(&num3,sizeof(float),1,fp); 

                    temp_surface.push_back(std::vector<Real>());

                    temp_surface[3*i+j].push_back( (Real)num1 );
                    temp_surface[3*i+j].push_back( (Real)num2 );
                    temp_surface[3*i+j].push_back( (Real)num3 );
                }

                fread(c,2,1,fp);
            }

            fclose(fp);
          
            makeWatertight(temp_surface);
        }
             
        // sort the generated pointlist, then merge vertices 
        void
        TriangulatedIF::makeWatertight( const std::vector<std::vector<Real> >& temp_surface) noexcept
        {
            std::vector< pointToElement > temp_ptlist;

            for(long i=0;i<this->normalList.size();++i)
            {
                faceList.push_back(Vec3ui(3*i, 3*i+1,3*i+2) );           

                for(int j=0;j<3;++j)
                {
                    temp_ptlist.push_back(pointToElement(temp_surface[3*i+j],i,j) );                        
                }
            }
            std::sort(temp_ptlist.begin(),temp_ptlist.end());

            mergeVertex(temp_ptlist);
        }
    
        // use a stupid algorithm here to merge vertices, acceleration may applied    
        void
        TriangulatedIF::mergeVertex(std::vector<pointToElement>& temp_ptlist) noexcept
        {
            long i=0,j=0,indexCounter=0,stepCounter;
        
            long mergeNumber=0;
        
            for(i=0;i<temp_ptlist.size();++i)
            {
                stepCounter=0;

                mergeNumber=0;
      
                vertList.push_back(Vec3r(temp_ptlist[i].physLocation[0],temp_ptlist[i].physLocation[1],temp_ptlist[i].physLocation[2])); 
            
                faceList[temp_ptlist[i].elementIndex].v[temp_ptlist[i].pointIndex] = indexCounter;
           
                for(j=i+1;j<temp_ptlist.size();++j)
                {
                    if(distance(temp_ptlist[i].physLocation,temp_ptlist[j].physLocation)<eps)
                    {
                        //merge vertex
                        faceList[temp_ptlist[j].elementIndex].v[temp_ptlist[j].pointIndex] = indexCounter;                    
                        temp_ptlist[j].mergeTag=1;                    
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
                    }
                }
                i+=stepCounter-1;
                if(stepCounter==0)
                {
                    break;
                }
            }
        }    
        
        void
        TriangulatedIF::buildDistance () noexcept
        {
            const Real* dx = geom().CellSize();
            const Real* plo = geom().ProbLo();
            for (MFIter mfi(distanceMF_); mfi.isValid(); ++mfi)
            {
                const Box& vbox = this->grids()[mfi.index()];
                Vec3r local_origin(plo[0] + vbox.smallEnd()[0]*dx[0],
                                   plo[1] + vbox.smallEnd()[1]*dx[1],
                                   plo[2] + vbox.smallEnd()[2]*dx[2]);
                Array3r phi_grid;
                double dx1 = double(dx[0]);
             
                make_level_set3(faceList, vertList, normalList,local_origin, dx1,
                                vbox.length(0),vbox.length(1),vbox.length(2), phi_grid);

                const auto& d = distanceMF_.array(mfi);
                const int* lo = vbox.loVect();
                const int* hi = vbox.hiVect();
                for (int k=lo[2]; k<=hi[2]; ++k)
                {
                    int kL=k-lo[2];
                    for (int j=lo[1]; j<=hi[1]; ++j)
                    {
                        int jL=j-lo[1];
                        for (int i=lo[0]; i<=hi[0]; ++i)
                        {
                            int iL=i-lo[0];
                            d(i,j,k) = phi_grid(iL,jL,kL);
                        }
                    }
                }
            }
            WriteSingleLevelPlotfile("Distance.out",distanceMF_,{"distance"},geom_,0.0,0);
        }
    }//endnamespace EB2
}//endnamespace amrex
