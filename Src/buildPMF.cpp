#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Print.H>

using namespace amrex;

const bool verbose_DEF = false;

static
void
PrintUsage (const char* progName)
{
    std::cout << '\n';
    std::cout << "Usage:" << '\n';
    std::cout << progName << '\n';
    std::cout << '\t' << "[premixfile = <filename>]" << '\n';
    std::cout << '\n';
    exit(1);
}


std::string
getFileRoot(const std::string& infile)
{
    std::vector<std::string> tokens = Tokenize(infile,std::string("."));
    std::string result = tokens[0];
    for (int i=1; i<tokens.size()-1; ++i)
        result = result + "." + tokens[i];
    return result;
}

void makeRoutine(std::vector<Vector<std::string> >& X,
                 std::string&                      routineName,
                 std::string&                      fileName);

int findIndex(const std::vector<std::string>& tokens, const std::string& str)
{
  int result = -1;
  for (int i=0; i<tokens.size() && result<0; ++i) {
    if (tokens[i] == str) result = i;
  }
  return result;
}

std::string ReplaceString(const std::string source_string,
                          const std::string old_substring,
                          const std::string new_substring)
{
    // Can't replace nothing.
    if (old_substring.empty())
        return source_string;

    // Find the first occurrence of the substring we want to replace.
    size_t substring_position = source_string.find(old_substring);

    // If not found, there is nothing to replace.
    if (substring_position == std::string::npos)
        return source_string;

    // Return the part of the source string until the first occurance of
    // the old substring + the new replacement substring + the result of the same function on the remainder.
    return source_string.substr(0,substring_position) +
      new_substring + ReplaceString(source_string.substr(substring_position+ old_substring.length(),
                                                         source_string.length()
                                                         - (substring_position + old_substring.length())),
                                    old_substring,
                                    new_substring);
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);
    {
    ParmParse pp;

    if (pp.contains("help"))
        PrintUsage(argv[0]);

    std::string premixfile = "";   pp.query("premixfile",premixfile);
    std::string chem1Dfile = "";   pp.query("chem1Dfile",chem1Dfile);
    std::string canterafile = "";   pp.query("canterafile",canterafile);

    bool use_premix = (premixfile != "");
    bool use_chem1D = (chem1Dfile != "");
    bool use_cantera = (canterafile != "");

    int filetypeCount=0;
    if (use_premix) filetypeCount++;
    if (use_chem1D) filetypeCount++;
    if (use_cantera) filetypeCount++;
    if (filetypeCount!=1)
      Abort("Must specific either a premix, chem1D or cantera file");
    
    std::string outfile="interp.f"; pp.query("outfile",outfile);
    std::string routineName="pmf"; pp.query("routineName",routineName);

    std::vector<Vector<std::string> > X;
    if (use_premix)
    {
      int Nspec=-1; pp.get("Nspec",Nspec);
      const int Nvals = Nspec + 4;
      
      const std::string TWOPNT("TWOPNT");
      const std::string FINAL("FINAL");

      std::cout << "Reading data from " << premixfile << "..." << std::endl;
      std::ifstream is(premixfile.c_str(),std::ios::in);
      if (!is.good())
        amrex::Abort(std::string((std::string("File ") + premixfile +
                                  std::string(" not opened"))).c_str());
      std::string line;
      std::getline(is,line);

      while (line.size()==0)
      {
        std::getline(is,line);
      }

      bool inData = false;
      bool finalDataSetRead = false;
      while (is.good() && !finalDataSetRead)
      {
        std::getline(is,line);

        std::vector<std::string> tokens;

        if (line.size() != 0)
          tokens = Tokenize(line,std::string(" :"));

        // In data, first value is integer pt number
        if (inData)
        {
          if (tokens.size() == Nvals + 1)
          {
            Vector<std::string> vals(Nvals);
            for (int i=0; i<Nvals; ++i)
              vals[i] = tokens[1+i];
            X.push_back(vals);
          }
          else
          {
            inData = false;
          }
        }

        if ((tokens.size() == Nvals))
        {
          X.clear();
          inData = true;
        }

        if ((tokens.size() == Nvals + 1) && (tokens[0] == std::string("X")))
        {
          amrex::Abort("In the data");
          X.clear();
          inData = true;
        }

        if (tokens.size()>=2 && tokens[0] == std::string("Total") && tokens[1] == std::string("CPUtime"))
          finalDataSetRead = true;
      }
      std::cout << "Read " << X.size() << " values" << std::endl;
    }
    else if (use_chem1D)
    {
      std::ifstream is(chem1Dfile.c_str(),std::ios::in);
      if (!is.good())
        amrex::Abort(std::string((std::string("File ") + premixfile +
                                  std::string(" not opened"))).c_str());

      int nPoints = -1;
      int nVars = -1;
      std::string line;

      // Find some integers in the header, the first in the list does not
      //  indicate an integer, but rather syas to stop processing the header
      Vector<std::string> str = {"FILE_STRUCTURE_COLUMNS_CONTAINING",
                                 "NUMBER_OF_GRIDPOINTS",
                                 "NUMBER_OF_VARIABLES",
                                 "NUMBER_OF_SPECIES"};
      std::map<std::string,int> res;

      int MassflowComp = -1;
      int Tcomp = -1;
      int Rcomp = -1;
      int Nspec = -1;

      {
        bool stopFound = false;
        std::getline(is,line);
        while (is.good() && !stopFound)
        {
          auto tokens = Tokenize(line,"[] ");
          if (tokens.size() > 0) {
            if (tokens[0] == str[0]) {
              stopFound = true;
            }
            else
            {
              for (int i=1; i<str.size(); ++i) {
                if (tokens[0] == str[i]) {
                  is >> res[str[i]];
                }
              }
            }
          }
          if (!stopFound) 
            std::getline(is,line);
        }

        nPoints = res["NUMBER_OF_GRIDPOINTS"];
        nVars = res["NUMBER_OF_VARIABLES"];
        Nspec = res["NUMBER_OF_SPECIES"];

        AMREX_ALWAYS_ASSERT(nPoints > 0);
        AMREX_ALWAYS_ASSERT(nVars > 0);
        AMREX_ALWAYS_ASSERT(Nspec > 0);
        
        std::getline(is,line);
        auto tokens = Tokenize(line," ");
        
        // Print() << "Species names: ";
        // for (int i=0; i<Nspec; ++i) {
        //   Print() << tokens[2+i] << " ";
        // }
        // Print() << std::endl;

        MassflowComp = -1;
        Tcomp = findIndex(tokens,"Temp");
        Rcomp = findIndex(tokens,"Density");
        MassflowComp  = findIndex(tokens,"Massflow");
      }

      X.resize(nPoints);
      for (int i=0; i<nPoints; ++i) {
        std::getline(is,line);

        line = ReplaceString(line,"E+","D+");
        line = ReplaceString(line,"E-","D-");
        
        auto tokens = Tokenize(line," ");
        AMREX_ALWAYS_ASSERT(tokens.size() == nVars + 2); // #, X, nVars

        // Set X vector from token data
        X[i].resize(Nspec + 4);
        X[i][0] = tokens[1]; // X
        X[i][1] = tokens[Tcomp];
        X[i][2] = std::to_string(std::stod(tokens[MassflowComp]) / std::stod(tokens[Rcomp]));
        X[i][3] = tokens[Rcomp];
        for (int j=0; j<Nspec; ++j) {
          X[i][4+j] = tokens[j+2];
        }
      }
    } else if (use_cantera) {
        // aja hack
        int Nspec=-1; pp.get("Nspec",Nspec);
	//
	// Read data
	//
        std::cout << "Reading " << canterafile << std::endl;
	std::ifstream is(canterafile.c_str());
        std::string line;
	std::vector<std::string> tokens;

	int nPts(0);

	// header
	std::getline(is,line);
	if (line.size() != 0) {
	  tokens = Tokenize(line,std::string(" ,"));
	  std::cout << "tokens = " << tokens.size() << std::endl;
	}

	do {
	  std::getline(is,line);
	  nPts++;
	} while (!is.eof());
	nPts--;  // empty last line

	std::cout << "nPts = " << nPts << std::endl;
	
	// destination
	// Nvals = 4 + nSpec: x / T / vel / rho / spec
	// stuff state into a fab
	X.resize(nPts);

	// source
	// 5 + nSpec: x / u / V / T / rho / X
	int nEntries=5+Nspec;
	int iZ=0;
	int iU=iZ+1;
	int iV=iU+1; // ???
	int iT=iV+1;
	int iR=iT+1;
	int iX=iR+1;
	Real vals[nEntries];

	// reset and read off header
	is.clear();
	is.seekg(0);
	std::getline(is,line);

	// loop over number of points and corresponding box
	for (int i=0; i<nPts; i++) {

	  // read line of data
	  std::getline(is,line);

	  // sanity check
	  if (line.size() == 0)
	    Abort("Empty line/EOF reading 1d data");
	  
	  // tokenize
	  tokens = Tokenize(line,std::string(","));

	  // if good, populate data
	  char tmpStr[32];
	  if (tokens.size()==nEntries) {
	    X[i].resize(Nspec + 4);
	    
	    for (int j=0; j<nEntries; j++) {
	      vals[j]=atof(tokens[j].c_str());
	    }
	    // space (m -> cm)
	    sprintf(tmpStr,"%.12e",vals[iZ] * 1.e2);
	    X[i][0] = tmpStr;
	    // temp
	    sprintf(tmpStr,"%.12e",vals[iT]);
	    X[i][1] = tmpStr;
	    // vel (m/s -> cm/s)
	    sprintf(tmpStr,"%.12e",vals[iU] * 1.e2);
	    X[i][2] = tmpStr;
	    // density (kg/m^3 -> g/cm^3)
	    sprintf(tmpStr,"%.12e",vals[iR] * 1.e-3);
	    X[i][3] = tmpStr;
	    // spec
	    for (int n=0; n<Nspec; n++) {
	      sprintf(tmpStr,"%.12e",vals[iX+n]);
	      X[i][4+n]=tmpStr;
	    }
	    
	  } else {
	    std::cout << "Expecting " << nEntries 
		      << " entries, but found " << tokens.size() << std::endl;
	    std::cout << line << std::endl;
	    Abort("Wrong number of entries reading 1d data");
	  }

	}

    }


    makeRoutine(X,routineName,outfile);
    }
    amrex::Finalize();
    return 0;
}

void makeRoutine(std::vector<Vector<std::string> >& X,
                 std::string&                      routineName,
                 std::string&                      fileName)
{
    std::ofstream os;
    os.open(fileName.c_str());

    int N=X.size();
    if (N==0)
    {
        amrex::Abort("No data to write");
    }

    int M=X.front().size(); // First one is position
    
    os << "c-------------------------------------------------\n";
    os << "c                                                 \n";
    os << "c                                                 \n";
    os << "c                              -USAGE-            \n";
    os << "c                                                 \n";
    os << "c     OVERVIEW: call " <<  routineName <<
        "(xlo,xhi,y,Niny)   \n";
    os << "c     computes y(xlo,xhi) by averaging data over  \n";
    os << "c     the range [xlo,xhi]. Here x is a scalar, y  \n";
    os << "c     a vector, with Niny components used. Piece- \n";
    os << "c     constant Extrapolation is applied to points \n";
    os << "c     laying outside the range of the orig data   \n";           
    os << "c                                                 \n";
    os << "c     domain: x in [" << X.front()[0] << "  \n";;
    os << "c                  ," << X.back()[0] << "]  \n";
    os << "c     number of states in vector y:" << M-1 << "        \n";
    os << "c     number of x data points:" << N << "             \n";
    os << "c                                                 \n";  
    os << "c     CALLING PROGRAM NEEDS: none                 \n";
    os << "c     PROCEDURE NEEDS:  none                      \n";
    os << "c     COMPILER NEEDS: f77, long names (14 char)   \n";
    os << "c     SYSTEM NEEDS: none                          \n";
    os << "c     NOTES: This routine was generated by        \n";
    os << "c     subroutine makeRoutine.                     \n";
    os << "c                                                 \n";
    os << "c-------------------------------------------------\n";
    os << "c                                                 \n";
    os << "c                          -SPECIFICATION-        \n";
    os << "c     SYNTAX:                                     \n";
    os << "c     subroutine " << routineName;
    os << "(xlo,xhi:dp,y:dp arr, Niny: int)\n";
    os << "c     returns: none                               \n";
    os << "c     REQUIRES: none                              \n";
    os << "c     MODIFIES: y, N                              \n";
    os << "c     SIGNALS: none                               \n";
    os << "c     EFFECTS: computes the average value of y    \n";
    os << "c     over the range [xlo,xhi] and returns number \n";
    os << "c     of components used, N.                      \n";
    os << "c                                                 \n";
    os << "c-------------------------------------------------\n";
    os << "c                                                 \n";
    os << "c                         -IMPLEMENTATION-        \n";
    os << "c     AUTHOR: auto-generated by makeAveInterp     \n";
    os << "c     MODIFICATION HISTORY: none                  \n";
    os << "c     NOTES: interpolation x,y data is defined in \n";
    os << "c     data statements for efficiency.             \n";
    os << "c     REQUIRES: none                              \n";
    os << "c     LOCAL VARIABLES:                            \n";
    os << "c     xlo,xhi.........(dp) average y on this range\n";
    os << "c     y_vector........(1D dp arr) computed y(x)   \n";
    os << "c     Niny............(int) no. of used components\n";
    os << "c     x_data..........(dp) x data points          \n";
    os << "c     y_data..........(dp) y(x) data points       \n";
    os << "c     N...............(int) number of x points    \n";
    os << "c     M...............(int) dimension of y_data   \n";
    os << "c     i,j,k...........(int) loop index            \n";
    os << "c     lo_{lo,hi}side..(int) idx of data near xlo  \n";          
    os << "c     hi_{lo,hi}side..(int) idx of data near xhi  \n";          
    os << "c     x1..............(dp) interpolation point    \n";         
    os << "c     x2..............(dp) interpolation point    \n";         
    os << "c     y1..............(dp) interpolation value    \n";         
    os << "c     y2..............(dp) interpolation value    \n";         
    os << "c     dydx............(dp) slope of linear interp \n";
    os << "c     ylo.............(dp) value at xlo           \n";         
    os << "c     yhi.............(dp) value at xhi           \n";         
    os << "c     GLOBAL VARIABLES:none                       \n";
    os << "c     PROCEDURES: none                            \n";
    os << "c     FILES: none                                 \n";
    os << "c     I/O: none                                   \n";
    os << "c-------------------------------------------------\n";
    os << "c                                                 \n";
    os << "c                             -CODE-              \n";
    os << "c     FIRST LINE:                                 \n";
    os << "      subroutine " << routineName;
    os << "(xlo,xhi,y_vector,Niny)\n";
      
    os << "c                                                 \n";
    os << "c     declare local variables                     \n";
    os << "c                                                 \n"; 
    os << "      integer N,M,i,j,k,lo_loside,lo_hiside       \n";
    os << "      integer hi_loside,hi_hiside,Niny            \n";
    os << "      double precision xlo,xhi,y_vector(*),sum    \n";
    os << "      double precision ylo,yhi,x1,y1,x2,y2,dydx   \n";
    os << "      double precision x_data(" << N << ")              \n";
    os << "      double precision y_data(" << N << "," << M-1 << ")       \n";
    os << "      data N,M /" << N << "," << M-1 << "/                   \n";

    // write number of components used in y_vector
    os << "      Niny = " << M-1 << "                                \n";
    // write interpolation data 
    os << "c                                                 \n";
    os << "c     interpolation x,y data                      \n";
    os << "c                                                 \n"; 

    for (int i=0; i<N; ++i)
    {
        const Vector<std::string>& data = X[i];
        os << "      data x_data(" << i+1 << ") /" <<  data[0]  << "/ \n";
        for (int j=1; j<M; ++j)
        {
            os << "      data y_data(" << i+1 << "," << j << ") /" << data[j]  << "/ \n";
        }
    }

    // write the interpolation routine
    os << "c                                                 \n";
    os << "c     interpolation routine                       \n";
    os << "c                                                 \n";
    os << "      lo_loside = 0\n";
    os << "      lo_hiside = 0\n";
    os << "      hi_loside = 0\n";
    os << "      hi_hiside = 0\n";
    os << "      if (xlo .le. x_data(1)) then\n";
    os << "         lo_loside = 1\n";
    os << "         lo_hiside = 1\n";
    os << "      end if\n";
    os << "      if (xhi .le. x_data(1)) then\n";
    os << "         hi_loside = 1\n";
    os << "         hi_hiside = 1\n";
    os << "      end if\n";
    os << "      if (xlo .ge. x_data(N)) then\n";
    os << "         lo_loside = N\n";
    os << "         lo_hiside = N\n";
    os << "      end if\n";
    os << "      if (xhi .ge. x_data(N)) then\n";
    os << "         hi_loside = N\n";
    os << "         hi_hiside = N\n";
    os << "      end if\n";
    os << "      if (lo_loside.eq.0) then\n";
    os << "         do i = 1, N-1                           \n";
    os << "            if ( (xlo .ge. x_data(i))\n";
    os << "     &           .and.\n";
    os << "     &           (xlo .le. x_data(i+1)) ) then\n";
    os << "               lo_loside  = i\n";
    os << "               lo_hiside  = i+1\n";
    os << "            end if\n";
    os << "         end do\n";
    os << "      end if\n";
    os << "      if (hi_loside.eq.0) then            \n";
    os << "         do i = 1, N-1                           \n";
    os << "            if ( (xhi .ge. x_data(i))\n";
    os << "     &           .and.\n";
    os << "     &           (xhi .le. x_data(i+1)) ) then\n";
    os << "               hi_loside = i\n";
    os << "               hi_hiside = i + 1\n";
    os << "            end if\n";
    os << "         end do\n";
    os << "      end if\n";
    os << "         \n";
    os << "      do j = 1, M\n";
    os << "         \n";
    os << "         x1 = x_data(lo_loside)\n";
    os << "         y1 = y_data(lo_loside,j)\n";
    os << "         \n";
    os << "         x2 = x_data(lo_hiside)\n";
    os << "         y2 = y_data(lo_hiside,j)\n";
    os << "          \n";
    os << "         if (lo_loside.eq.lo_hiside) then\n";
    os << "            dydx = 0.d0\n";
    os << "         else\n";
    os << "            dydx = (y2-y1)/(x2-x1)\n";
    os << "         end if\n";
    os << "         \n";
    os << "         ylo = y1 + dydx*(xlo - x1)\n";
    os << "         \n";
    os << "         if (lo_loside .eq. hi_loside) then\n";
    os << "            \n";
    os << "            yhi = y1 + dydx*(xhi - x1)\n";
    os << "            \n";
    os << "            y_vector(j) = 0.5d0*(ylo + yhi)\n";
    os << "            \n";
    os << "         else\n";
    os << "            \n";
    os << "            sum = (x2 - xlo) * 0.5d0 * (ylo + y2)\n";
    os << "            \n";
    os << "            x1 = x_data(hi_loside)\n";
    os << "            y1 = y_data(hi_loside,j)\n";
    os << "         \n";
    os << "            x2 = x_data(hi_hiside)\n";
    os << "            y2 = y_data(hi_hiside,j)\n";
    os << "            \n";
    os << "            if (hi_loside.eq.hi_hiside) then\n";
    os << "               dydx = 0.d0\n";
    os << "            else\n";
    os << "               dydx = (y2-y1)/(x2-x1)\n";
    os << "            end if\n";
    os << "            \n";
    os << "            yhi = y1 + dydx*(xhi - x1)\n";
    os << "            \n";
    os << "            sum = sum + (xhi - x1)*0.5d0*(yhi+y1)\n";
    os << "            \n";
    os << "            do k = lo_hiside,hi_loside-1\n";
    os << "               \n";
    os << "               sum = sum + (x_data(k+1)-x_data(k))\n";
    os << "     &              * 0.5d0\n";
    os << "     &              * (y_data(k,j) + y_data(k+1,j))\n";
    os << "               \n";
    os << "            end do\n";
    os << "            \n";
    os << "            y_vector(j) = sum / (xhi - xlo)\n";
    os << "            \n";
    os << "         end if\n";
    os << "      end do\n";
    os << "      end                                         \n";
    os << "c     LAST LINE                                   \n"; 
    os << "c-------------------------------------------------\n";
}

