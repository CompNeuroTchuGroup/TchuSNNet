#ifndef RandomConnectivity_hpp
#define RandomConnectivity_hpp

class Synapse;
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>

#include <fstream>

class RandomConnectivity : public Connectivity {
protected:

    NeuronInt noSourceNeurons{}; //no of inputs per neuron (synapses for which each neuron is postsynaptic)
    virtual void SetNoSourceNeurons(NeuronInt readNoSourceNeurons);

public:

    RandomConnectivity(Synapse* synapse ,const GlobalSimInfo*  infoGlobal);
    ~RandomConnectivity() override = default;

    virtual void        ConnectNeurons() override;
    double       GetExpectedConnections() const override {return noSourceNeurons;}
    std::string   GetTypeStr() const override {return IDstringRandomConnectivity;}
    double              GetConnectionProbability() const;

    virtual void SaveParameters(std::ofstream& wParameterStream, std::string idString) const override;
    virtual void LoadParameters(const std::vector<FileEntry>& parameters) override;
};

#endif /* RandomConnectivity_hpp */
