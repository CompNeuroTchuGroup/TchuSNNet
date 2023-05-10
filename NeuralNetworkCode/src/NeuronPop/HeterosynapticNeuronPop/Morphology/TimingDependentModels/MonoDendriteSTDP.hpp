//
// Created by Saif Ahmed on 27.06.21.
//

#ifndef NEURALNETWORK_MONODENDRITE_H
#define NEURALNETWORK_MONODENDRITE_H

#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "../Morphology.hpp"
#include "../SynapseSpines/CoopSynapseSpine.hpp"

typedef std::shared_ptr<BaseSynapseSpine> BaseSpinePtr;

class HeteroCurrentSynapse;

class MonoDendriteSTDP: public Morphology {
protected:

    //Moved from Morphology.hpp
    std::vector<bool> spikedSynapses;
    std::vector<std::pair<unsigned long, double>> preSpikes;
    std::vector<double> postSpikes;
    std::vector<std::pair<unsigned long, double>> theta_changes;
    std::vector<std::pair<unsigned long, double>> weight_changes;
    //End
    double tauTheta{}; // decay constant of heterosynaptic effects in spines
    double lambdaDist{}; // decay constant of heterosynaptic effects over distance between synapses
    double tauDelay{}; // decay constant of heterosynaptic effects over inter-synapse spike timing difference

    double thetaExpDecay{1.0};
    double dendriticLength{}; // this would change in case of more complex dendritic geometry (atm it is a single 1D dendrite)
    double synapticGap{}; //  minimum gap between syanpses along dendrite

    double preFactorLTP{};
    double preFactorLTD{};

    long synapseIdGenerator{}; // variable used to allocate new synapses
    double nextPos{};
    bool stepWeights{};
    std::vector<unsigned long> weightStepBoundary{};
    std::vector<double> weightStepValue{};
    unsigned long currWightStepId{};

    std::vector<bool> integratePostSpike{};
    std::vector<bool> integratePreSpike{};

    std::vector<std::shared_ptr<CoopSynapseSpine>> synapseDataCoop;

    double initialWeights{1.0};

    double base_ltp{};
    double base_ltd{};

    double incr_ltp{};
    double decr_ltd{};

    void thetaDecay();
    void updateCooperativity(unsigned long synId, unsigned long neighborId);
    void pseudoCoop(unsigned long synId, unsigned long neighborId);

    virtual void updateLTP(unsigned long synId) = 0;
    virtual void updateLTD(unsigned long synId) = 0;

    virtual double gLTP(double deltaT) const = 0;
    virtual double gLTD(double deltaT) const = 0;

    virtual double aLTP(double theta) const = 0;
    virtual double aLTD(double theta) const = 0;

    virtual double getDistanceEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const = 0;
    virtual double getTimingEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const = 0;

public:
    explicit MonoDendriteSTDP(GlobalSimInfo* info);
    ~MonoDendriteSTDP() override = default;

    void advect() override;
    void Reset() override;
    void RecordPostSpike() override;
    void RecordExcitatoryPreSpike(int spikedSynapseId) override;

    void SaveParameters(std::ofstream * stream, std::string neuronPreId) override;
    void LoadParameters(std::vector<std::string> *input) override;

    virtual BaseSpinePtr AllocateNewSynapse(const HeteroCurrentSynapse& synapse) override;


    //Revirtualization

};


#endif //NEURALNETWORK_MONODENDRITE_H
