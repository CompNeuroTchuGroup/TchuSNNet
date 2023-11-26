//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef NEURALNETWORK_MORPHOLOGY_H
#define NEURALNETWORK_MORPHOLOGY_H

#include "../Morphology/SynapseSpines/BaseSynapseSpine.hpp"
#include "./../../GlobalFunctions.hpp"
#include <cmath>
#include <limits>
#include <vector>

#include <algorithm>
#include <iostream>
#include <string>

enum WeightNormalization {
    NOPNormalization,
    HardNormalization,
    SoftMaxNormalization
};

class Morphology {

  protected:


    GlobalSimInfo *infoGlobal;

    std::vector<BaseSpinePtr> baseSpineData;


    // double weightsSum {};
    double totalPostSpikes{};
    double totalPreSpikes{};
    bool   postSpiked{false};

    bool         distributeWeights{false};
    bool         userSeed{false};
    int          seed{0};
    std::mt19937 generator;


    double              minWeight{0.0};
    double              initialWeights{1.0};
    double              maxWeight{2.0};

    virtual void Reset() = 0;


  public:
    explicit Morphology(GlobalSimInfo *infoGlobal);
    virtual ~Morphology() = default;

    virtual void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const;
    virtual void LoadParameters(const std::vector<FileEntry> &parameters);
    virtual void CheckParameters(const std::vector<FileEntry> &parameters);

    virtual BaseSpinePtr AllocateNewSynapse(BranchTargeting &synapse) = 0;
    double GenerateSynapticWeight(); // Here we generate the synaptic weight to be allocated when a synapse is allocated
    int    GetNoSynapses() const { return static_cast<int>(baseSpineData.size()); }
    virtual std::string GetType() const = 0;

    virtual void Advect() = 0;
    virtual void RecordPostSpike();
    virtual void RecordExcitatoryPreSpike(BaseSpinePtr spinePtr);
    // Getters
    std::vector<double>         GetIndividualSynapticProfile(signed long synapseId) const;
    std::string                 GetIndividualSynapticProfileHeaderInfo() const;
    virtual std::vector<double> GetOverallSynapticProfile() const;
    virtual std::string         GetOverallSynapticProfileHeaderInfo() const;
    // virtual void CalcMorphoPlasticityEvents() {return;};
    // friend std::vector<signed long> getSpikedSynapsesFromMorphology(const Morphology&); // This function is not
    // necessary as the spikedSynapses is not used outside of the class
    //  signed long GetSynapseCount() const;
    virtual int    GetMaxGapDelay(int delayPerMicroMeter) { return 0; }
    virtual double GetSynapticDistanceToSoma(int synapseId) = 0;
    virtual void   PostConnectSetUp(){};

    virtual bool IgnoreJDParameters() const { return false; }
};

#endif // NEURALNETWORK_MORPHOLOGY_H
