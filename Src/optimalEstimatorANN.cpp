#include <optimalEstimatorANN.H>

#include <cstring>
#include <fstream>
#include <vector>

Net::Net(int input_size, amrex::Vector<int> neurons, int output_size)
{
  n_layers = neurons.size();
  layers.resize(n_layers + 1, nullptr);
  if (n_layers == 0) {
    layers[0] =
      register_module("fc", torch::nn::Linear(input_size, output_size));
    layers[0]->to(torch::kDouble);
  } else {
    layers[0] =
      register_module("fc0", torch::nn::Linear(input_size, neurons[0]));
    layers[0]->to(torch::kDouble);
    layers[n_layers] = register_module(
      "fc" + std::to_string(n_layers),
      torch::nn::Linear(neurons[n_layers - 1], output_size));
    layers[n_layers]->to(torch::kDouble);
    for (int n = 1; n < n_layers; n++) {
      layers[n] = register_module(
        "fc" + std::to_string(n),
        torch::nn::Linear(neurons[n - 1], neurons[n]));
      layers[n]->to(torch::kDouble);
    }
  }
  initializeWeights();
}

// Implement the forward pass
torch::Tensor
Net::forward(torch::Tensor x)
{
  // following Berger et al. (2018), linear activation function, after tansig
  // function
  for (int n = 0; n < n_layers; n++) {
    x = activate(layers[n]->forward(x));
  }
  x = layers[n_layers]->forward(x);

  return x;
}

torch::Tensor
Net::activate(torch::Tensor x)
{
  return torch::tanh(x);
}

// Function to initialize weights
void
Net::initializeWeights()
{
  // Calculate the scaling factor
  for (int n = 0; n <= n_layers; n++) {
    torch::nn::init::xavier_uniform_(layers[n]->weight);
  }
}

// ---------------------------------------------------------------------------
// Normalisation bounds file; see optimalEstimatorANN.H for the layout.
// ---------------------------------------------------------------------------

namespace {

constexpr char minmaxMagic[8] = {'P', 'e', 'l', 'e', 'O', 'E', 'M', 'M'};
constexpr std::int32_t minmaxVersion = 1;

void
writeName(std::ostream& os, const std::string& s)
{
  const std::int32_t n = static_cast<std::int32_t>(s.size());
  os.write(reinterpret_cast<const char*>(&n), sizeof(n));
  os.write(s.data(), n);
}

std::string
readName(std::istream& is, const std::string& path)
{
  std::int32_t n = 0;
  is.read(reinterpret_cast<char*>(&n), sizeof(n));
  // A corrupt length would otherwise ask for an arbitrary allocation.
  if (!is || n < 0 || n > 1024) {
    amrex::Abort("Malformed variable name in " + path);
  }
  std::string s(static_cast<std::size_t>(n), '\0');
  is.read(s.data(), n);
  if (!is) {
    amrex::Abort("Truncated variable name in " + path);
  }
  return s;
}

std::string
join(const amrex::Vector<std::string>& v)
{
  std::string s;
  for (std::size_t i = 0; i < v.size(); i++) {
    s += (i ? " " : "") + v[i];
  }
  return s;
}

std::string
join(const amrex::Vector<int>& v)
{
  std::string s;
  for (std::size_t i = 0; i < v.size(); i++) {
    s += (i ? " " : "") + std::to_string(v[i]);
  }
  return s;
}

} // namespace

void
writeMinMax(
  const std::string& path,
  const amrex::Vector<std::string>& features,
  const amrex::Vector<std::string>& targets,
  const amrex::Vector<int>& neurons,
  const amrex::Vector<amrex::Real>& f_min,
  const amrex::Vector<amrex::Real>& f_max,
  const amrex::Vector<amrex::Real>& t_min,
  const amrex::Vector<amrex::Real>& t_max)
{
  std::ofstream file(path, std::ios::binary);
  if (!file.good()) {
    amrex::Abort("Could not open minmax file " + path + " for writing");
  }

  const std::int32_t ver = minmaxVersion;
  const std::int32_t nF = static_cast<std::int32_t>(features.size());
  const std::int32_t nT = static_cast<std::int32_t>(targets.size());
  file.write(minmaxMagic, sizeof(minmaxMagic));
  file.write(reinterpret_cast<const char*>(&ver), sizeof(ver));
  file.write(reinterpret_cast<const char*>(&nF), sizeof(nF));
  file.write(reinterpret_cast<const char*>(&nT), sizeof(nT));
  for (const auto& s : features) {
    writeName(file, s);
  }
  for (const auto& s : targets) {
    writeName(file, s);
  }

  const std::int32_t nL = static_cast<std::int32_t>(neurons.size());
  file.write(reinterpret_cast<const char*>(&nL), sizeof(nL));
  for (const int n : neurons) {
    const std::int32_t w = static_cast<std::int32_t>(n);
    file.write(reinterpret_cast<const char*>(&w), sizeof(w));
  }

  // Always double, so that a FLOAT build can read a DOUBLE build's file.
  auto put = [&](const amrex::Vector<amrex::Real>& v) {
    const std::vector<double> d(v.begin(), v.end());
    file.write(
      reinterpret_cast<const char*>(d.data()), sizeof(double) * d.size());
  };
  put(f_min);
  put(f_max);
  put(t_min);
  put(t_max);

  file.close();
  if (!file) {
    amrex::Abort("Failed to write minmax file " + path);
  }
}

void
readMinMax(
  const std::string& path,
  const amrex::Vector<std::string>& features,
  const amrex::Vector<std::string>& targets,
  const amrex::Vector<int>& neurons,
  amrex::Vector<amrex::Real>& f_min,
  amrex::Vector<amrex::Real>& f_max,
  amrex::Vector<amrex::Real>& t_min,
  amrex::Vector<amrex::Real>& t_max)
{
  const int nF = static_cast<int>(features.size());
  const int nT = static_cast<int>(targets.size());
  f_min.resize(nF);
  f_max.resize(nF);
  t_min.resize(nT);
  t_max.resize(nT);

  std::ifstream file(path, std::ios::binary);
  if (!file.good()) {
    amrex::Abort("Could not open minmax file " + path);
  }

  char magic[sizeof(minmaxMagic)] = {};
  file.read(magic, sizeof(magic));
  const bool hasHeader =
    file && std::memcmp(magic, minmaxMagic, sizeof(magic)) == 0;

  if (!hasHeader) {
    // Written before the header existed. Nothing in such a file can be
    // checked, so the old failure modes still apply to it.
    file.clear();
    file.seekg(0);
    amrex::Print()
      << "\n*** WARNING: " << path
      << " carries no header and was written by an older\n"
      << "    optimalEstimatorTraining. The variable names it was built from "
         "and the\n"
      << "    precision it was written in cannot be checked, so a mismatch "
         "with this run\n"
      << "    would go unnoticed. Retrain to get a self-describing file.\n\n";

    auto getRaw = [&](amrex::Vector<amrex::Real>& v) {
      file.read(
        reinterpret_cast<char*>(v.dataPtr()), sizeof(amrex::Real) * v.size());
    };
    getRaw(f_min);
    getRaw(f_max);
    getRaw(t_min);
    getRaw(t_max);
    if (!file) {
      amrex::Abort(
        "Short read from " + path +
        ": it does not match the requested number of features and targets");
    }
    return;
  }

  std::int32_t ver = 0, nFfile = 0, nTfile = 0;
  file.read(reinterpret_cast<char*>(&ver), sizeof(ver));
  file.read(reinterpret_cast<char*>(&nFfile), sizeof(nFfile));
  file.read(reinterpret_cast<char*>(&nTfile), sizeof(nTfile));
  if (!file) {
    amrex::Abort("Truncated header in minmax file " + path);
  }
  if (ver != minmaxVersion) {
    amrex::Abort(
      path + " is format version " + std::to_string(ver) +
      ", but this build writes and reads version " +
      std::to_string(minmaxVersion) + ". Retrain to regenerate it.");
  }
  if (nFfile != nF || nTfile != nT) {
    amrex::Abort(
      path + " was written for " + std::to_string(nFfile) + " feature(s) and " +
      std::to_string(nTfile) + " target(s), but " + std::to_string(nF) +
      " feature(s) and " + std::to_string(nT) +
      " target(s) were requested. Training and inference must be given the "
      "same lists.");
  }

  amrex::Vector<std::string> fileFeatures(nF), fileTargets(nT);
  for (int n = 0; n < nF; n++) {
    fileFeatures[n] = readName(file, path);
  }
  for (int n = 0; n < nT; n++) {
    fileTargets[n] = readName(file, path);
  }
  // Order matters as much as membership: the bounds are applied per component.
  if (fileFeatures != features) {
    amrex::Abort(
      path + " was trained with features '" + join(fileFeatures) + "', but '" +
      join(features) +
      "' were requested. Training and inference must be given the same list, "
      "in the same order.");
  }
  if (fileTargets != targets) {
    amrex::Abort(
      path + " was trained with targets '" + join(fileTargets) + "', but '" +
      join(targets) +
      "' were requested. Training and inference must be given the same list, "
      "in the same order.");
  }

  std::int32_t nL = 0;
  file.read(reinterpret_cast<char*>(&nL), sizeof(nL));
  if (!file || nL < 0 || nL > 1024) {
    amrex::Abort("Malformed layer count in " + path);
  }
  amrex::Vector<int> fileNeurons(nL);
  for (std::int32_t n = 0; n < nL; n++) {
    std::int32_t w = 0;
    file.read(reinterpret_cast<char*>(&w), sizeof(w));
    fileNeurons[n] = static_cast<int>(w);
  }
  if (!file) {
    amrex::Abort("Truncated layer list in " + path);
  }
  // torch::load does not object to a checkpoint whose layers do not match the
  // network it is loaded into, so this is the only place the architecture is
  // checked at all.
  if (fileNeurons != neurons) {
    amrex::Abort(
      path + " was trained with neurons '" + join(fileNeurons) + "', but '" +
      join(neurons) +
      "' were requested. The architecture is rebuilt from this argument before "
      "the checkpoint is loaded and must match the training run.");
  }

  auto get = [&](amrex::Vector<amrex::Real>& v) {
    std::vector<double> d(v.size());
    file.read(reinterpret_cast<char*>(d.data()), sizeof(double) * d.size());
    for (std::size_t i = 0; i < d.size(); i++) {
      v[i] = static_cast<amrex::Real>(d[i]);
    }
  };
  get(f_min);
  get(f_max);
  get(t_min);
  get(t_max);
  if (!file) {
    amrex::Abort("Short read from " + path + ": the file is truncated.");
  }
}
