#include <analysis_util.H>

#include <AMReX_ParallelDescriptor.H>
#include <fstream>

namespace analysis_util {

std::string
parse_title(std::istream& is)
{
  std::string line;
  std::getline(is, line);
  return line;
}

std::vector<std::string>
parse_var_names(std::istream& is)
{
  std::string line;
  std::getline(is, line);
  return amrex::Tokenize(line, std::string(", "));
}

MEFData
read_mef(const std::string& infile)
{
  std::ifstream ifs(infile.c_str());
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
    ifs.good(), "analysis_util::read_mef: cannot open file " + infile);

  MEFData data;
  data.title = parse_title(ifs);
  data.var_names = parse_var_names(ifs);

  ifs >> data.n_elts >> data.nodes_per_elt;

  data.nodes.readFrom(ifs);

  const int n_conn = data.n_elts * data.nodes_per_elt;
  data.connectivity.resize(n_conn);
  ifs.read(
    reinterpret_cast<char*>(data.connectivity.data()), sizeof(int) * n_conn);

  return data;
}

void
write_mef(const std::string& outfile, const MEFData& data)
{
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::ofstream ofs(
      outfile.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      ofs.good(), "analysis_util::write_mef: cannot open file " + outfile);

    ofs << data.title << "\n";
    for (int i = 0; i < static_cast<int>(data.var_names.size()); ++i) {
      if (i > 0)
        ofs << " ";
      ofs << data.var_names[i];
    }
    ofs << "\n";

    ofs << data.n_elts << " " << data.nodes_per_elt << "\n";

    data.nodes.writeOn(ofs);

    ofs.write(
      reinterpret_cast<const char*>(data.connectivity.data()),
      sizeof(int) * data.connectivity.size());
  }
  // All ranks wait until IOProcessor has finished writing so that
  // subsequent read_mef calls on all ranks see the complete file.
  amrex::ParallelDescriptor::Barrier();
}

} // namespace analysis_util
