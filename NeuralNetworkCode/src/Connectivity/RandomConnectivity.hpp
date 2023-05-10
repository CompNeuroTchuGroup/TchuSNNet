#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>
#include <valarray>
#include <fstream>
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"

#ifndef RandomConnectivity_hpp
#define RandomConnectivity_hpp

class Synapse;

class RandomConnectivity : public Connectivity
{
protected:

    unsigned long noSourceNeurons; //no of inputs per neuron (synapses for which each neuron is postsynaptic)
    virtual void SetNoSourceNeurons(unsigned long noSN);

public:

    RandomConnectivity(Synapse *syn,GlobalSimInfo  * info);
    ~RandomConnectivity()=default;

    virtual void        ConnectNeurons() override;
    unsigned long       GetNumberAverageSourceNeurons() override {return noSourceNeurons;}
    const std::string   GetTypeStr() override {return str_randomConnectivity;}
    double              GetConnectionProbability();

    virtual void SaveParameters(std::ofstream * stream, std::string id_str) override;
    virtual void LoadParameters(std::vector<std::string> *input) override;
};

#endif /* RandomConnectivity_hpp */
