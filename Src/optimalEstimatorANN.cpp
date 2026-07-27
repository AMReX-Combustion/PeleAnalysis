#include <optimalEstimatorANN.H>

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
