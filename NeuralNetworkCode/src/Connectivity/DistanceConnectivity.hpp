#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>
#include <valarray>
#include <fstream>
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"

#ifndef DistanceConnectivity_hpp
#define DistanceConnectivity_hpp

class Synapse;

class DistanceConnectivity : public Connectivity
{
protected:

    double PeakProbability;
	double StdProbability;
	double a;
	int exact;
public:

	DistanceConnectivity(Synapse *syn,GlobalSimInfo  * info);
    ~DistanceConnectivity()=default;

    void                ConnectNeurons();
	void                ConnectNeuronsExact();
	unsigned long		GetNumberAverageSourceNeurons();
    const std::string         GetTypeStr() override {return str_distanceConnectivity;}

    void SaveParameters(std::ofstream * stream,std::string id_str);
    void LoadParameters(std::vector<std::string> *input);
};

#endif /* DistanceConnectivity_hpp */
