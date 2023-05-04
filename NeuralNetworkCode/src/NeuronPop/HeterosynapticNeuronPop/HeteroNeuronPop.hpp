//
// Created by Saif Ahmed on 03.08.21.
//

#ifndef NEURALNETWORK_HETEROPOP_H
#define NEURALNETWORK_HETEROPOP_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>

#include "Morphology/TimingDependentModels/MonoDendriteSTDPTazerart.hpp"
#include "Morphology/TimingDependentModels/MonoDendriteSTDPTazerartRelative.hpp"
#include "Morphology/TimingDependentModels/MonoDendriteSTDPBiWindow.hpp"
#include "Morphology/BranchedResourceSTDPModel/BranchedResourceHeteroSTDP.hpp"
#include "../../GlobalFunctions.hpp"
#include "../NeuronPop.hpp"

typedef std::shared_ptr<BaseSynapseSpine> BaseSpinePtr;

class HeteroCurrentSynapse;
class SynapseSpine;
class HeteroNeuronPop: public NeuronPop {

protected:
    std::vector<std::unique_ptr<Morphology>> morphology{};
    bool isBranched{false}; //PARAMETERLOADING with Set isBranchedTrue

public:
    explicit HeteroNeuronPop(GlobalSimInfo * info, int id);
    ~HeteroNeuronPop() override = default;

    void SaveParameters(std::ofstream * stream) override;
    void LoadParameters(std::vector<std::string> *input) override;

    void PostConnectSetUp() override;
    //Getters
    std::valarray<double> GetIndividualSynapticProfile(unsigned long neuronId, unsigned long synapseId)  override;
    std::valarray<double> GetOverallSynapticProfile(unsigned long neuronId) override;
    std::string GetIndividualSynapticProfileHeaderInfo() const override;
    std::string GetOverallSynapticProfileHeaderInfo() const override;
    unsigned long GetSynapseCount(unsigned long neuronId);
    double GetWeight(unsigned long neuronId, unsigned long synapseId);

    //Setters
    void SetBranchedTrue(){isBranched=true;}

    // Heterosynaptic functionality
    BaseSpinePtr AllocateNewSynapse(unsigned long neuronId, HeteroCurrentSynapse& syn) override;
    void RecordExcitatorySynapticSpike(unsigned long neuronId, unsigned long synapseId) override;


    //Inheritance optimizations
    bool HasHeterosynapticPlasticity() override {return true;}
    bool IsBranchedBool() override {return isBranched;}

    unsigned long GetNumberOfSynapses() const override {return morphology.size();}
    // Testing, check the purpose of this
    //friend std::vector<unsigned long> getSpikedSynapses(const HeteroNeuronPop&, unsigned long neuronId); //currently not in use

    //Morphology& getMorphologyRef(int neuronId); //Return a reference to Morpho object

    // STDP Analysis
    //void triggerStatOut(std::string dirPath); //currently not in use
    //void printThetasAndWeights();  //currently not in use



};


#endif //NEURALNETWORK_HETEROPOP_H
