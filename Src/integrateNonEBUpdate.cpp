/*
 * =====================================================================================
 *
 *       Filename:  integrateNonEBUpdate
 *
 *    Description:  Integral tool that takes variables from a plotfile and returns volume integral of that variable(s)-can run in parallel
 *
 *        Version:  1.0
 *        Created:  03/10/2020 01:37:30 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Victor H. Zendejas Lopez 
 *   Organization:  CCSE @LBL
 *
 * =====================================================================================
 */


#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

using namespace amrex;
using namespace std;

static void PrintUsage (char* progName)
{
   cout << "\nUsage:\n"
         << progName
         << "\n\tinfile = inputFileName comps=[n1 n2 ...]"
         << "\n\t[-help]"
         << "\n\n";
    exit(1);
}


static Real SumThisComp(AmrData &amrData, string whichVar , int finestLevel)
{

   Real sum(0.0);

   // Get levels of refinement 
   int nLevels = finestLevel+1;
   // Intersect variablenumber (must come after stoichiometry)
   Vector<unique_ptr< MultiFab>> mf(nLevels);
   for (int iLevel=0; iLevel < nLevels; iLevel++){
      Print()<< "Reading Level " << iLevel << "\n";
      int ngrow(0);

      // Make space for each component and one for the intersect flag
      DistributionMapping dm(amrData.boxArray(iLevel));
      mf[iLevel].reset(new MultiFab(amrData.boxArray(iLevel), dm, 1, ngrow));
      amrData.FillVar(*mf[iLevel], iLevel, whichVar);      
   }

   for (int iLevel=0; iLevel<finestLevel; iLevel++) {
      BoxArray baf = mf[iLevel+1]->boxArray();
      baf.coarsen(amrData.RefRatio()[iLevel]);   
 
      for (MFIter mfi(*mf[iLevel]); mfi.isValid(); ++mfi) {
         FArrayBox& myFab = (*mf[iLevel])[mfi]; 
         int idx = mfi.index();
         const Box& bx = mfi.validbox();

         vector < pair < int,Box > > isects = baf.intersections(mf[iLevel]->boxArray()[idx]); 
         for (int ii = 0; ii < isects.size(); ii++){
            myFab.setVal(0,isects[ii].second,0,1);
         }

         sum +=myFab.sum(0,1);
      }
   }

   ParallelDescriptor::ReduceRealSum(sum);
   amrData.FlushGrids(amrData.StateNumber(whichVar));
   return sum;
}


int main(int argc, char* argv[])
{
    Initialize(argc,argv);
    {
    if(argc == 1){
    //   PrintUsage(argv[0]);
    }     
   
    ParallelDescriptor::StartParallel(&argc, &argv);
    ParmParse pp;

    string Output_Variables;
    //Count the number of terms in the Output_Varibales 
    int nOutput_Variables = (pp.countval("Output_Variables"));
    Vector < string > Out_Var(nOutput_Variables); // Creates a vector Out_var length  
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
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);    
    DataServices plotDataType(Plot_File, fileType);
     
    if (!plotDataType.AmrDataOk()) // If nothing is stored within plotDataType then run this loop
       DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = plotDataType.AmrDataRef(); 
    cout << "   ...done." << endl;
 
    // Most likely, the variable was averaged down conservatively, lev 0 may be all that is reqd
    int finestLevel = amrData.FinestLevel();

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

     if (ParallelDescriptor::IOProcessor()){
       cout << "TIME = " << setw(15) << amrData.Time();
    }

    for (int v = 0; v<whichVar.size(); ++v){
          Real sumComp = SumThisComp(amrData, whichVar[v] ,finestLevel);
          if (ParallelDescriptor::IOProcessor()){
             cout << " int(" << whichVar[v] << ") = "
              << sumComp << '\n';
       }
    }
      }
       Finalize();
       return 0;
}
