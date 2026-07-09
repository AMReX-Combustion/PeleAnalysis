#ifndef analysis_util_H
#define analysis_util_H

#include <istream>
#include <string>
#include <vector>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

namespace analysis_util {

// --- AMR utilities ---
[[nodiscard]]
amrex::Vector<amrex::iMultiFab>
get_covered_mf(
  const amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<int>&             ref_ratios);

// --- Integration tools ---
[[nodiscard]]
std::unique_ptr<amrex::Gpu::ManagedVector<amrex::Real>>
integrate(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  const amrex::Vector<amrex::Geometry>& geoms,
  const amrex::Vector<int>&             ref_ratios,
  int                                   scomp,
  int                                   ncomp,
  const amrex::Vector<int>&             axes_to_integrate);

// --- Gradient tools ---
void gradient(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  const amrex::Vector<amrex::Geometry>& geoms,
  int                                   scomp,
  int                                   ncomp,
  amrex::Vector<amrex::MultiFab>&       grad_mf);

// --- String / file utilities ---

// Strip directory prefix: "path/to/plt00001" -> "plt00001"
std::string get_file_root(const std::string& infile);

// Read a single title/label line from a stream
std::string parse_title(std::istream& is);

// Read and tokenize a variable names line from a stream
std::vector<std::string> parse_var_names(std::istream& is);

// Find component index of a named variable; aborts if not found and
// abort_if_not_found=true
int find_var_index(
  const amrex::Vector<std::string>& var_names,
  const std::string&                name,
  bool                              abort_if_not_found = true);

// --- I/O tools ---

// Read the `n_files` ParmParse option (if present) and apply it to AMReX so
// that subsequent plotfile writes use at most `n_files` data files
// (VisMF::SetNOutFiles). Useful to limit file counts for large post-processing
// campaigns. Returns the number of output files now in effect. When `n_files`
// is not specified, the current AMReX default is left unchanged.
int set_plot_nfiles();

struct PlotfileData
{
  amrex::Vector<amrex::MultiFab> mf;
  amrex::Vector<amrex::Geometry> geoms;
  amrex::Vector<int>             ref_ratios;
  amrex::Real                    time;
  int                            n_lev;
  amrex::Vector<std::string>     var_names; // all variables present in the file
};

// Load selected variables from an AMReX plotfile into MultiFabs + Geometries.
[[nodiscard]]
PlotfileData read_plotfile(
  const std::string&                infile,
  const amrex::Vector<std::string>& var_names,
  int                               finest_level = 1000,
  int                               n_grow = 0,
  const amrex::Vector<int>&         is_per = {});

// Write multi-level MultiFabs to an AMReX plotfile.
// ref_ratios defaults to all-2 if empty.
void write_plotfile(
  const std::string&                    outfile,
  const amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<std::string>&     var_names,
  const amrex::Vector<amrex::Geometry>& geoms,
  amrex::Real                           time = 0.0,
  const amrex::Vector<int>&             ref_ratios = {});

// --- MEF file tools ---

struct MEFData
{
  std::string              title;
  std::vector<std::string> var_names;
  amrex::FArrayBox         nodes;        // component-major layout [nNodes x nComp]
  amrex::Vector<int>       connectivity; // flat array [n_elts * nodes_per_elt]
  int                      n_elts;
  int                      nodes_per_elt;
};

[[nodiscard]]
MEFData read_mef(const std::string& infile);

void write_mef(const std::string& outfile, const MEFData& data);

// --- PlotfileData convenience overloads ---

inline amrex::Vector<amrex::iMultiFab>
get_covered_mf(const PlotfileData& d)
{
  return get_covered_mf(d.mf, d.ref_ratios);
}

inline std::unique_ptr<amrex::Gpu::ManagedVector<amrex::Real>>
integrate(
  const PlotfileData&       d,
  int                       scomp,
  int                       ncomp,
  const amrex::Vector<int>& axes)
{
  return integrate(d.mf, d.geoms, d.ref_ratios, scomp, ncomp, axes);
}

inline void
gradient(
  const PlotfileData&            d,
  int                            scomp,
  int                            ncomp,
  amrex::Vector<amrex::MultiFab>& grad_mf)
{
  return gradient(d.mf, d.geoms, scomp, ncomp, grad_mf);
}

inline void
write_plotfile(const std::string& outfile, const PlotfileData& d)
{
  write_plotfile(outfile, d.mf, d.var_names, d.geoms, d.time, d.ref_ratios);
}

} // namespace analysis_util

#endif // analysis_util_H
