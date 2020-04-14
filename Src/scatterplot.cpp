/*
 * =====================================================================================
 *
 *       Filename:  scatterplot.cpp
 *
 *    Description:  Create an output file from a Plot file that can be used to generate plots
 *
 *        Version:  1.0
 *        Created:  02/12/2020 02:51:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Victor H. Zendejas Lopez
 *   Organization:  CCSE
 *
 * =====================================================================================
 */


/*This code works in series only. Running in Parrelel will cause
Errors to accumulate in outputfiles */

#include <iostream>
#include <string>
#include <fstream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

using namespace amrex;
using namespace std;


int main (int argc, char* argv[])
{
   
   Initialize(argc,argv);
   {

      // ParmParse is a class that has stored functions  
      ParmParse pp; 
      // declares Output_Variables as an string-comes from the input file
      string Output_Variables;
      //Count the number of terms in the Output_Varibales 
      int nOutput_Variables = (pp.countval("Output_Variables"));
      Vector < string > Out_Var(nOutput_Variables); // Creates a vector Out_var length nOU..es 
      // Loop over Output_Variables to extract information from each element in Out_Var

      for(int ivarOut = 0; ivarOut < nOutput_Variables; ++ivarOut){ 
         pp.get("Output_Variables", Out_Var[ivarOut],ivarOut);
      }
      string Output_Format; 
      pp.get("Output_Format", Output_Format); 
      string Plot_File; 
      pp.get("Plot_File", Plot_File);

      string Output_File; 
      pp.get("Output_File", Output_File); 
      // Extract Plot data
      // 1. DataService-Class-setBatchMode is function in that class
      // 2. FileType is a class, plotType is a variable that takes in arguments Amrvis::NEWPLT
      //3. DataServices creates a variable that takesin Plot_File and plotType as an arguement 
      //4. amrData is a variable that contains the directory of where the data from 
      //the plot files are located. AmrDataRef is a fucntion acting on plotDataType.
    
      DataServices::SetBatchMode();
      Amrvis::FileType plotType(Amrvis::NEWPLT); 
      DataServices plotDataType(Plot_File, plotType);

      if (!plotDataType.AmrDataOk()) // If nothing is stored within plotDataType then run this loop
         DataServices::Dispatch(DataServices::ExitRequest, NULL);

      AmrData& amrData = plotDataType.AmrDataRef(); 
      cout << "   ...done." << endl;

      Vector<string> whichVar(nOutput_Variables);
      // Read in variable list - adjusting for stoichiometry

      for(int v = 0; v < nOutput_Variables; v++) {
         pp.get("Output_Variables", whichVar[v], v);
      }
      // Check the names of the variables are present in the plotfile
      Vector < int > destFills(nOutput_Variables);

      for (int v = 0; v < nOutput_Variables; v++){
         destFills[v] = v;

         if (amrData.StateNumber(whichVar[v])<0) {
            string message="Bad variable name ("+whichVar[v]+")";
            Error(message.c_str());
         }
      } 
      // Get domain size
      Vector<Real> probLo=amrData.ProbLo();
      Vector<Real> probHi=amrData.ProbHi();
      // Get probDomain
      Vector<Box> probDomain = amrData.ProbDomain();
      // Get levels of refinement 
      int finestLevel = amrData.FinestLevel();    	
      int nLevels = finestLevel+1;
      // Intersect variable number (must come after stoichiometry)
      int isVar=nOutput_Variables;

      // Load the data       
      // 
      // 1. boxArray is a class for storing collection of boxes on a single AMR level
      // 2. DistributionMapping is a class that describes which process owns the data
      // on the domain specified by the boxes in boxArray
      //
      //  Vector<unique_ptr<MultiFab>> mf(nLevels);
      Vector<unique_ptr< MultiFab>> mf(nLevels);
      for (int iLevel=0; iLevel < nLevels; iLevel++){
         Print()<< "Reading Level " << iLevel << "\n";
         int ngrow(0);

         // Make space for each component and one for the intersect flag
         DistributionMapping dm(amrData.boxArray(iLevel));
         mf[iLevel].reset(new MultiFab(amrData.boxArray(iLevel), dm, nOutput_Variables+1, ngrow));
         // Put the value 1 in the intersect flag
         mf[iLevel]->setVal(1.,isVar,1);
         // Load the data (make a copy of whichVar to make sure it's the same size as destFills) 

         Vector<string> loadWhichVar(nOutput_Variables);
         for (int v = 0; v < nOutput_Variables; v++){        
            loadWhichVar[v] = whichVar[v];
         } 
         amrData.FillVar(*mf[iLevel], iLevel, loadWhichVar, destFills);      
      }

      // Set the intersect flag to zero if necessary
      for (int iLevel=0; iLevel<finestLevel; iLevel++) {
         BoxArray baf = mf[iLevel+1]->boxArray();
         baf.coarsen(amrData.RefRatio()[iLevel]);   
         Print()<< "Tagging Level  " << iLevel << "\n";
        
         for (MFIter mfi(*mf[iLevel]); mfi.isValid(); ++mfi) {
            FArrayBox& myFab = (*mf[iLevel])[mfi]; 
            int idx = mfi.index();
            vector < pair<int,Box> > isects = baf.intersections(mf[iLevel]->boxArray()[idx]);
            for (int ii = 0; ii < isects.size(); ii++){
               myFab.setVal(0,isects[ii].second,isVar,1);
            }
         }
      }
//      for (int iLevel=0; iLevel < nLevels; iLevel++){
//         VisMF::Write(*mf[iLevel],"MFvis_Level"+to_string(iLevel));
//      }

      // Create an outputfile to write data to 
      ofstream ofs;
      ofs.open(Output_File);
     
     // Load Header file for output file 
      Vector<string> loadWhichVar(nOutput_Variables);
      for (int v = 0; v < nOutput_Variables; v++){        
         loadWhichVar[v] = whichVar[v];
         ofs << "   " << loadWhichVar [v];
      }
      ofs << "\n";

      for (int iLevel=0; iLevel < nLevels; iLevel++){

         for (MFIter mfi(*mf[iLevel]); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            const auto& Tagfab = mf[iLevel]->array(mfi,isVar);
            const auto& DataFab = mf[iLevel]->array(mfi);
            
            // Write data to output file
            
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(box,i,j,k,
                  {

                  if( Tagfab(i,j,k) > 0.5){
                  for(int n = 0; n < nOutput_Variables; ++n){
                     ofs << " " <<  DataFab(i,j,k,n);
                  }
                  ofs << "\n";
                  }
                  });
         } 
      }
      ofs.close();
   }
   Finalize();
   return 0;
}

