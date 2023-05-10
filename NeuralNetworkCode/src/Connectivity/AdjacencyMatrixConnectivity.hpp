#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>
#include <valarray>
#include <fstream>
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"

#ifndef AdjacencyMatrixConnectivity_hpp
#define AdjacencyMatrixConnectivity_hpp

class Synapse;

class AdjacencyMatrixConnectivity : public Connectivity
{
protected:

    unsigned long noSourceNeurons; //no of inputs per neuron (synapses for which each neuron is postsynaptic)
    std::vector<std::vector<int>> connectivity_matrix;
    std::string connectionProbFile;
    void SetNoSourceNeurons(unsigned long noSN);

public:

    AdjacencyMatrixConnectivity(Synapse* syn, GlobalSimInfo* info);
    ~AdjacencyMatrixConnectivity()=default;

    void                ConnectNeurons();
    unsigned long       GetNumberAverageSourceNeurons() { return noSourceNeurons; }
    const std::string         GetTypeStr() override { return str_adjacencyMatrixConnectivity; }
    double              GetConnectionProbability();
    void                GetConnectionWeightsFromFile(std::string filepath);

    void SaveParameters(std::ofstream* stream, std::string id_str);
    void LoadParameters(std::vector<std::string>* input);

};

#endif /* AdjacencyMatrixConnectivity_hpp */
