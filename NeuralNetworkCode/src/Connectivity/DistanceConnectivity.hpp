#ifndef DistanceConnectivity_hpp
#define DistanceConnectivity_hpp

#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"
#include "../Synapse/Synapse.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>

#include <fstream>

class Synapse;

class DistanceConnectivity : public Connectivity
{
protected:

    double peakProbability{};
	double stdProbability{1};
	double projectionLengthFactor{1};//==a
	int isConnectionExact{};
public:

	DistanceConnectivity(Synapse* synapse,const GlobalSimInfo*  infoGlobal);
    ~DistanceConnectivity() override = default;

    void                ConnectNeurons();
	void                ConnectNeuronsExact();
	double		GetExpectedConnections() const override;
    std::string         GetTypeStr() const override {return IDstringDistanceConnectivity;}

    void SaveParameters(std::ofstream& wParameterStream,std::string idString) const;
    void LoadParameters(const std::vector<FileEntry>& parameters);
};

#endif /* DistanceConnectivity_hpp */
