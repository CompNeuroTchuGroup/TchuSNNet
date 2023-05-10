#ifndef NEURALNETWORK_MORPHOLOGY_H
#define NEURALNETWORK_MORPHOLOGY_H


#include "./SynapseSpines/BaseSynapseSpine.hpp"
#include "./../../../GlobalFunctions.hpp"
#include "../../../Synapse/HeteroCurrentSynapse.hpp"
#include <vector>
#include <cmath>
#include <limits>
#include <valarray>
#include <algorithm>
#include <iostream>
#include <string>

typedef std::shared_ptr<BaseSynapseSpine> BaseSpinePtr;
typedef std::shared_ptr<Branch> BranchPtr;

enum WeightNormalization {
    NOPNormalization, HardNormalization, SoftMaxNormalization
};

class Morphology {

protected:
    void hardNormalize();
    void softMaxNormalize();

    // STDP Analysis
    std::vector<std::pair<unsigned long, double>> preSpikes;
    std::vector<double> postSpikes;

    std::vector<std::pair<unsigned long, double>> theta_changes;
//    std::vector<double> theta_changes;
//    std::vector<double> theta_changes;

    std::vector<std::pair<unsigned long, double>> weight_changes;

    GlobalSimInfo * info;

    std::vector<BaseSpinePtr> baseSynapseData;

    double weightsSum {};
    double totalPostSpikes {};
    double totalPreSpikes {};
    bool postSpiked{false};

    bool distributeWeights{false};
    int seed{0};
    std::mt19937 generator{ std::mt19937(seed) };

    double lastPostSpikeTime;
    std::vector<bool> spikedSynapses;
    std::vector<int> spikedSynapsesId;

    WeightNormalization weightNormalization {NOPNormalization};
    double minWeight {0.0};
    double initialWeights {1.0};
    double maxWeight {2.0};

    bool decayWeights {false};
    double WeightDecayConstant{1.0};
    double weightExpDecay {};

    unsigned long totalPlasticityEvents{};

    virtual void Reset()=0;
    virtual void NormalizeWeights();

    virtual void WeightDecay();

public:
    explicit Morphology(GlobalSimInfo * info);
    virtual ~Morphology() = default;

    virtual void SaveParameters(std::ofstream * stream, std::string neuronPreId);
    virtual void LoadParameters(std::vector<std::string> *input);

    virtual BaseSpinePtr AllocateNewSynapse(const HeteroCurrentSynapse& synapse)=0;
    double GenerateSynapticWeight();// Here we generate the synaptic weight to be allocated when a synapse is allocated
    virtual const std::string GetType() = 0;

    virtual void advect() = 0;
    virtual void RecordPostSpike();
    virtual void RecordExcitatoryPreSpike(int spikedSynapseId);
    //Getters
    std::valarray<double> GetIndividualSynapticProfile(unsigned long synapseId) const;
    std::string GetIndividualSynapticProfileHeaderInfo() const;
    virtual std::valarray<double> GetOverallSynapticProfile();
    virtual std::string GetOverallSynapticProfileHeaderInfo() const;
    //virtual void CalcMorphoPlasticityEvents() {return;};
    //friend std::vector<unsigned long> getSpikedSynapsesFromMorphology(const Morphology&); // This function is not necessary as the spikedSynapses is not used outside of the class
    unsigned long GetSynapseCount() const;
    double GetWeight(unsigned long synapseId) const;
    virtual void PostConnectSetUp(){};
    // STDP Analysis
    //void triggerStatOut(std::string dirPath);

    //void printThetasAndWeights();

    //Qualifying methods
    virtual bool const IsBranchedBool() {return false;}

};

//std::vector<unsigned long> getSpikedSynapsesFromMorphology(const Morphology&);

#endif //NEURALNETWORK_MORPHOLOGY_H
