#ifndef Connectivity_hpp
#define Connectivity_hpp

#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>

#include <fstream>
#include "../GlobalFunctions.hpp"

class Synapse;

class Connectivity{
protected:

    const GlobalSimInfo* infoGlobal;
    bool                        userSeed{false};//by default, seed can be changed
    int                         seed{};
    std::mt19937  generator;

    Synapse*  synapse;
    //synapticTargets has to be changed to synapticTargets!!!!
    //std::vector<std::vector<signed long>> synapticTargets; //the list with postsynaptic (or target) neurons for each neuron of the presynaptic population {array of pointer to vectors}

public:

    Connectivity(Synapse* synapse,const GlobalSimInfo*  infoGlobal);
    virtual ~Connectivity() = default;

    // TODO/Suggestion: Replace pointer return with a return by constant reference. YES PLEASE
    // https://github.com/saiftyfirst/BP_Demos/blob/master/C%2B%2B/constRef_vs_pointer.cpp
    //const std::vector<signed long>& GetTargetList(long preNeuronId) const {return synapticTargets.at(preNeuronId);}
    virtual double GetExpectedConnections() const = 0;
    virtual std::string GetTypeStr() const = 0;

    virtual void ConnectNeurons() = 0;

    // virtual void WriteConnectivity(const std::string& filename,signed long noNeuronsConnectivity);
    // virtual void WriteDistributionD(const std::string& filename, signed long noNeuronsDelay);
    // virtual void WriteDistributionJ(const std::string& filename, signed long noNeuronsJPot);

    virtual void SaveParameters(std::ofstream& wParameterStream, std::string idString) const;
    virtual void LoadParameters(const std::vector<FileEntry>& parameters);


    void SetSeed(std::mt19937& seedGenerator);
    // void SetSeed(int readSeed);

    // TODO/Suggestion: Replace pointer with raw value. For primitives there is no advantage returning a pointer unless it is changed by the caller

    // virtual void SetDistributionD();

    // TODO/Suggestion: Replace pointer with raw value. For primitives there is no advantage returning a pointer unless it is changed by the caller

    // virtual void SetDistributionJ();

    // void Test();

    //Move to synapse, including types to avoid calls
    // const std::vector<std::pair<signed long, BaseSpinePtr>>& GetSynapticTargets(const signed long sourceNeuron) const;
    // virtual double GetDistributionJ(long preNeuronId, long postNeuronId) const;
    // virtual int GetDistributionD(long preNeuronId, long postNeuronId);
    // int GetSeed() const {return seed;}
    // void SortSynapticTargetList();
};

#endif /* Connectivity_hpp */
