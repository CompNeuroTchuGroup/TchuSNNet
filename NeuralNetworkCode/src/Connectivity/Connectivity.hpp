#ifndef Connectivity_hpp
#define Connectivity_hpp

#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>
#include <valarray>
#include <fstream>
#include "../GlobalFunctions.hpp"

class Synapse;

class Connectivity{
protected:

    GlobalSimInfo               * info;
    bool                        fixSeed{false};
    int                         seed{};
    std::mt19937  generator;

    Synapse                     * synapse;
    std::vector<unsigned long>  * target_id; //the list with postsynaptic (or target) neurons for each neuron of the presynaptic population {array of pointer to vectors}
    std::vector<int>            * D_distribution; // the list with delays D that are associated with the synapses betweeen the source and its target neurons
    std::vector<double>         * J_distribution; // list of synaptic weights J
//	bool						HasPot;
	bool						HasDdistribution{false};
    bool                        HasJDistribution{false};
public:

    Connectivity(Synapse *syn,GlobalSimInfo  * info);
    virtual ~Connectivity();

    // TODO/Suggestion: Replace pointer return with a return by constant reference. YES PLEASE
    // https://github.com/saiftyfirst/BP_Demos/blob/master/C%2B%2B/constRef_vs_pointer.cpp
    std::vector<unsigned long>* GetTargetList(long preNeuronId){return &target_id[preNeuronId];}
    virtual unsigned long GetNumberAverageSourceNeurons() =0;
    virtual const std::string GetTypeStr() = 0;

    virtual void ConnectNeurons() = 0;
    virtual void WriteConnectivity(const std::string& filename,unsigned long noNeuronsConnectivity);
    virtual void WriteDistributionD(const std::string& filename, unsigned long noNeuronsDelay);
    virtual void WriteDistributionJ(const std::string& filename, unsigned long noNeuronsJPot);

    virtual void SaveParameters(std::ofstream * stream, std::string id_str);
    virtual void LoadParameters(std::vector<std::string> *input);

    int GetSeed() const { return seed; }
    void SetSeed(std::mt19937 *seedGenerator);
    void SetSeed(int s);

    // TODO/Suggestion: Replace pointer with raw value. For primitives there is no advantage returning a pointer unless it is changed by the caller
    virtual int* GetDistributionD(long preNeuronId, long postNeuronId);
    virtual void SetDistributionD();

    // TODO/Suggestion: Replace pointer with raw value. For primitives there is no advantage returning a pointer unless it is changed by the caller
    virtual double* GetDistributionJ(long preNeuronId, long postNeuronId);
    virtual void SetDistributionJ();

    void Test();

    //Removing dynamic_casts from main loop
    virtual std::vector<std::pair<unsigned long, unsigned long>>& getSynapticTargets(const unsigned long&);

};

#endif /* Connectivity_hpp */
