//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _MORPHOLOGY_BASE_CLASS_HPP
#define _MORPHOLOGY_BASE_CLASS_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "../Morphology/SynapseSpines/BaseSynapseSpine.hpp"
#include "./../../GlobalFunctions.hpp"

enum WeightNormalization { NOPNormalization, HardNormalization, SoftMaxNormalization };

class Morphology {
  protected:
    GlobalSimInfo *infoGlobal;

    std::vector<BaseSpinePtr> baseSpineData;

    // double weightsSum {};
    double totalPostSpikes {};
    double totalPreSpikes {};
    bool   postSpiked { false };

    bool         distributeWeights { false };
    bool         userSeed { false };
    int          seed { 0 };
    std::mt19937 generator;

    double minWeight { 0.0 };
    double initialWeights { 1.0 };
    double maxWeight { 2.0 };

    virtual void Reset() = 0;

    std::mutex _presynSpikeMutex;

  public:
    explicit Morphology(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters);
    virtual ~Morphology() = default;

    virtual void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const;
    virtual void LoadParameters(const std::vector<FileEntry> &parameters);
    virtual void CheckParameters(const std::vector<FileEntry> &parameters);

    virtual BaseSpinePtr AllocateNewSynapse(BranchTargeting &synapse) = 0;
    double               GenerateSynapticWeight();  // Here we generate the synaptic weight to be allocated when a synapse is allocated
    SynInt               GetNoSynapses() const { return static_cast<SynInt>(baseSpineData.size()); }
    virtual std::string  GetType() const = 0;

    virtual void Advect() = 0;
    virtual void RecordPostSpike();
    virtual void RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) = 0;
    // Record functions
    std::vector<double> GetIndividualSynapticProfile(signed long synapseId) const;
    std::string         GetIndividualSynapticProfileHeaderInfo() const;

    virtual std::vector<double> GetOverallSynapticProfile() const;
    virtual std::string         GetOverallSynapticProfileHeaderInfo() const;

    virtual std::vector<std::pair<std::string, double>> GetSteadyStateData() const;
    // virtual void CalcMorphoPlasticityEvents() {return;};
    // friend std::vector<signed long> getSpikedSynapsesFromMorphology(const Morphology&); // This function is not
    // necessary as the spikedSynapses is not used outside of the class
    //  signed long GetSynapseCount() const;
    // virtual int    GetMaxGapDelay(int delayPerMicroMeter) { return 0; }
    // virtual double GetSynapticDistanceToSoma(int synapseId) = 0;
    virtual void PostConnectSetUp() {};

    virtual bool IgnoreJcoupling() const { return false; }
    virtual bool HasSteadyState() const { return false; }
};

#endif  // NEURALNETWORK_MORPHOLOGY_H
