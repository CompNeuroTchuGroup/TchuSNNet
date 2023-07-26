//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_
#define _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_

#include "BranchedSynapseSpine.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <unordered_map>

class ResourceSynapseSpine : public BranchedSynapseSpine {

    protected:
    //Resource model N and K.
    // double kBasal{1.0};
    // double nBasal{1.0};

    // double kStimulus{0.0};
    // double nStimulus{0.0};


    // double kStimulusExpDecay{1.0};//This is of stimulus variable
    // double nStimulusExpDecay{1.0};//This is of stimulus variable

    //Resource model alpha
    double alphaBasal{1.0};
    double alphaStim{0.0}; //This is delta alpha
    double alphaStimBump{1.0};
    double alphaStimulusExpDecay{1.0};

    //Central vars
    double alphaResources{};

    //Morpho-copies
    double potDepRatio{1.0};


// If I have to go back here, an unordered_map<index, pair<count<pair<double, double>>>>, and find how to delete or vector
    // std::vector<double> kStimulusTempEffect{};//These need to be here for non-matrix management. 
    // std::vector<double> nStimulusTempEffect{};
    // std::vector<int> stimulusEffectCount{};
    // std::list<std::pair<double, int>> potAlphaTempCount;
    // std::vector<double> depAlphaTempAndCount;

    // int maxCount{100};

    // bool depressionFlagSTDP{false};
    // bool potentiationFlagSTDP{false};
    // bool updated{false};
    
    public:
    ResourceSynapseSpine() = default;
    ~ResourceSynapseSpine() override = default;
    //Getters
    double GetAlphaResources() const {return alphaResources;}
    // bool GetDepressionFlagSTDP(){return depressionFlagSTDP;}
    // bool GetPotentiationFlagSTDP(){return potentiationFlagSTDP;}
    // bool GetUpdatedFlag(){return updated;}
    //int GetMaxCount(){return maxCount;}//Not necessary for now
    //Setters
    void SetAlphaBasal(double alphaBasalInput){this->alphaBasal=alphaBasalInput;}
    void SetAlphaStim(double alphaStimInput){this->alphaStim=alphaStimInput;}
    //void SetAlphaStimulusRest(double alphaStimulusRest){alphaStimulus=alphaStimulusRest;}
    void SetAlphaExpDecay(double alphaStimExpDecay){this->alphaStimulusExpDecay=alphaStimExpDecay;}
    void SetAlphaStimBump(double alphaStimBump){this->alphaStimBump=alphaStimBump;}
    // void SetNBasal(int nBasalInput){nBasal=nBasalInput;}
    // void SetKBasal(int kBasalInput){kBasal=kBasalInput;}
    // void SetKExponentialDecay(int kStimulusExpDecayCalc){kStimulusExpDecay=kStimulusExpDecayCalc;}//Input is supposed to be exp(-dt/tau)
    // void SetNExponentialDecay(int nStimulusExpDecayCalc){nStimulusExpDecay=nStimulusExpDecayCalc;}
    // void SetMaxCount(int maxCountnew) {maxCount=maxCountnew;} //UNRESOLVED  Call in AllocateNewSynapse
    void SetPotentiationRatio(double ratio) {this->potDepRatio=ratio;}
    // void SetDepressionFlag(bool booleanFlag) {depressionFlagSTDP=booleanFlag;}
    // void SetPotentiationFlag(bool booleanFlag) {potentiationFlagSTDP=booleanFlag;}
    // void SetUpdatedFlag(bool flag){updated=flag;}
    //Alpha methods
    void ComputeAlphaResources();//Should be recalced at least once in the Spine Setup
    void DecayAlphaResources();
    void ComputeWeight(double resourceFactor);
    //Temp effets methods
    // void AddTempResourcesToSpine(double alphaStimulusInput);
    // bool ApplyAllTemporaryEffectsOnPostspike(const std::unordered_map<int, double>& STDPdecayMap);//input must be -1 if depression or STDPratio if potentiation, or the inverse swapping everything
    //void ApplyAllTempEffectsOnConflictPotentiation(double PotentiationDepressionRatio);//input must be -1 if depression or STDPratio if potentiation, or the inverse swapping everything
    // bool ApplyAllTemporaryEffectsOnDepression(double expDecayFactor);//input must be -1 if depression or STDPratio if potentiation, or the inverse swapping everything
    //Trace methods
    void AddTraceToAlpha(double computedTrace);
    //Stimulus vector methods
    // void TickStimulusCounts();//Called in Reset(), but both should be mutually exclusive with AATE above
    // void CullStimulusVectors();//Called in Reset()
    //Profile methods
    std::vector<double> GetIndividualSynapticProfile() const override;
    std::string GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif