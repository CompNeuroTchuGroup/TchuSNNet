#ifndef PoissonConnectivity_hpp
#define PoissonConnectivity_hpp
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>

#include <fstream>


class Synapse;

class PoissonConnectivity : public Connectivity {
protected:
    signed long  totalConnections{};
    double connectionLambda{};

public:

    PoissonConnectivity(Synapse* synapse,const GlobalSimInfo*  infoGlobal);
    ~PoissonConnectivity() override = default;

    void ConnectNeurons() override;
    std::string GetTypeStr() const override { return IDstringPoissonConnectivity; }

    double GetExpectedConnections() const override;

    void SaveParameters(std::ofstream& wParameterStream, std::string idString) const override;
    void LoadParameters(const std::vector<FileEntry>& parameters) override;
};

#endif /* PoissonConnectivity_hpp */
