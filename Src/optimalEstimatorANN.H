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

#endif
