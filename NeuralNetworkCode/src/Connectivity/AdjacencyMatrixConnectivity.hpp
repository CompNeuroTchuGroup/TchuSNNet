#ifndef _ADJACENCY_MATRIX_CONNECTIVITY_HPP
#define _ADJACENCY_MATRIX_CONNECTIVITY_HPP
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"
#include <iostream>
#include <random>
#include <vector>

class Synapse;

class AdjacencyMatrixConnectivity : public Connectivity {
protected:
  NeuronInt                     noSourceNeurons{}; // no of inputs per neuron (synapses for which each neuron is postsynaptic)
  std::vector<std::vector<int>> connectivityMatrix;
  std::string                   connectionProbFile;
  void                          SetNoSourceNeurons(NeuronInt readNoSourceNeurons);

public:
  AdjacencyMatrixConnectivity(Synapse *synapse, GlobalSimInfo *infoGlobal);
  ~AdjacencyMatrixConnectivity() override = default;

  void        ConnectNeurons();
  double      GetExpectedConnections() const override { return noSourceNeurons; }
  std::string GetTypeStr() const override { return IDstringAdjacencyMatrixConnectivity; }
  double      GetConnectionProbability() const;
  void        GetConnectionWeightsFromFile(std::string filepath);

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const;
  void LoadParameters(const std::vector<FileEntry> &parameters);
};

#endif /* AdjacencyMatrixConnectivity_hpp */
