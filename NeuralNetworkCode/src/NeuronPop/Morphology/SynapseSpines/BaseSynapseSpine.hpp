//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _BASE_SYNAPSE_SPINE_CLASS_HEADER_
#define _BASE_SYNAPSE_SPINE_CLASS_HEADER_

#include "../../../GlobalFunctions.hpp"

#include <string>
#include <memory>
#include <mutex>

class BaseSynapseSpine {

    protected:
    std::mutex _relativeCouplingMutex;
    //Legacy variables
    signed long idInMorpho{}; // id for synapse within its population
    NeuronInt preNeuronId{};
    NeuronInt postNeuronId{};
    PopInt prePopId{};
    double couplingStrength{1};
    double weight{}; //The negative weight comes from J, weight is just a factor to multiply (for now)
    //double lastSpike{};

    //Branched variables
    bool isBranchedBool{false};// For now it is not useful
    //Synapse related


    //Necessary mutex

    public:
    //Constructors
    BaseSynapseSpine() = default;
    virtual ~BaseSynapseSpine() = default;
    //BaseSynapseSpine(double weight, double lastSpike);// Not currently in use
    //Methods
    //Getters
    NeuronInt GetPreNeuronId() const {return preNeuronId;};
    NeuronInt GetPostNeuronId() const {return postNeuronId;};
    double GetWeight() const {return weight*couplingStrength;};
    double GetWeightUncoupled() const {return weight;};
    signed long GetIdInMorpho() const {return idInMorpho;};
    bool GetBranchedBool() const { return isBranchedBool;}
    //Setters
    void SetPreNeuronId(signed long neuronId){preNeuronId=neuronId;};
    void SetPostNeuronId(signed long neuronId){postNeuronId=neuronId;};
    void SetWeight(double weightIn){weight=weightIn;};
    void SetIdInMorpho(signed long idIn){idInMorpho=idIn;};
    void SetRelativeCouplingStrength(double rCouplingStrength);
    void SetPreNeuronPop(PopInt population){prePopId=population;}
    //Misc
    void AddToWeight(double change){weight+=change;}

    //Recorder functions
    virtual std::vector<double> GetIndividualSynapticProfile() const = 0;
    virtual std::string GetIndividualSynapticProfileHeaderInfo() const = 0;
};

#endif