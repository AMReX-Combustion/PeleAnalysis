#ifndef ANN_H
#define ANN_H

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <torch/torch.h>
#include <tuple>
#include <random>
#include <string>
#include <iostream>

// optimalEstimatorANN.cpp

struct Net : torch::nn::Module
{
public:
  // Constructor declaration
  Net(int input_size, amrex::Vector<int> neurons, int output_size);

  // Forward pass declaration
  torch::Tensor forward(torch::Tensor x);

private:
  // Forward pass declaration
  // torch::Tensor tansig(torch::Tensor x);
  torch::Tensor activate(torch::Tensor x);

  // Function to initialize weights
  void initializeWeights();

  // Declare layers
  // torch::nn::Linearfc1{nullptr}, fc2{nullptr}, fc3{nullptr}; // Add more
  // layers if needed

  int n_layers;
  amrex::Vector<torch::nn::Linear> layers;
};

// Normalisation bounds file, "<minmax_path>.bin".
//
// Written by optimalEstimatorTraining and read back by optimalEstimatorInfer,
// which has to de-normalise the network output with exactly the bounds the
// network was trained under. The payload is always double, like the network
// checkpoint, so the two executables may be built with different PRECISION.
//
// The file opens with a header naming the format version, the feature and
// target names and the hidden-layer widths, and readMinMax checks all of them
// against what the caller asked for. Getting the argument lists out of step
// between the two runs is the easy mistake to make - the tools take them on the
// command line and share nothing - and none of it is caught anywhere else:
// torch::load accepts a checkpoint whose layers do not match the network it is
// loaded into, so a wrong "neurons" runs to completion and writes a
// plausible-looking plotfile full of nonsense.
//
// Layout:
//   char[8]  magic "PeleOEMM"
//   int32    format version
//   int32    number of features, int32 number of targets
//   per name, features then targets: int32 length, then that many chars
//   int32    number of hidden layers, then that many int32 widths
//   double   f_min[nF], f_max[nF], t_min[nT], t_max[nT]
//
// Files written before the header existed are still read, with a warning; they
// are a bare dump of the four arrays in the writing build's own precision.
void writeMinMax(
  const std::string& path,
  const amrex::Vector<std::string>& features,
  const amrex::Vector<std::string>& targets,
  const amrex::Vector<int>& neurons,
  const amrex::Vector<amrex::Real>& f_min,
  const amrex::Vector<amrex::Real>& f_max,
  const amrex::Vector<amrex::Real>& t_min,
  const amrex::Vector<amrex::Real>& t_max);

void readMinMax(
  const std::string& path,
  const amrex::Vector<std::string>& features,
  const amrex::Vector<std::string>& targets,
  const amrex::Vector<int>& neurons,
  amrex::Vector<amrex::Real>& f_min,
  amrex::Vector<amrex::Real>& f_max,
  amrex::Vector<amrex::Real>& t_min,
  amrex::Vector<amrex::Real>& t_max);

#endif
