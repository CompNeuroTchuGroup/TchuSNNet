#ifndef _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_
#define _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_

#include "BranchedSynapseSpine.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <unordered_map>

typedef std::unordered_map<int, double> DHashMap;
typedef std::unordered_map<int, DHashMap> NestedDHashMap;
typedef std::pair<double, int> PairDI;

class ResourceSynapseSpine : public BranchedSynapseSpine {

    protected:
    //Resource model N and K.
    // double kBasal{1.0};
    // double nBasal{1.0};

    // double kStimmulus{0.0};
    // double nStimmulus{0.0};


    // double kStimmulusExpDecay{1.0};//This is of stimmulus variable
    // double nStimmulusExpDecay{1.0};//This is of stimmulus variable

    //Resource model alpha
    double alphaBasal{1.0};
    double alphaStimmulus{0.0}; //This is delta alpha
    double alphaStimmulusExpDecay{1.0};

    //Central vars
    double alphaResources{};

    //Morpho-copies
    double PotentiationDepressionRatio{1.0};

// If I have to go back here, an unordered_map<index, pair<count<pair<double, double>>>>, and find how to delete or vector
    // std::vector<double> kStimmulusTempEffect{};//These need to be here for non-matrix management. 
    // std::vector<double> nStimmulusTempEffect{};
    // std::vector<int> stimmulusEffectCount{};
    std::list<PairDI> potentiationAlphaTempAndCount;
    std::vector<double> depressionAlphaTempAndCount;

    int maxCount{100};

    // bool depressionFlagSTDP{false};
    // bool potentiationFlagSTDP{false};
    bool updated{false};
    
    public:
    //Getters
    double GetAlphaResources(){return alphaResources;}
    // bool GetDepressionFlagSTDP(){return depressionFlagSTDP;}
    // bool GetPotentiationFlagSTDP(){return potentiationFlagSTDP;}
    bool GetUpdatedFlag(){return updated;}
    //int GetMaxCount(){return maxCount;}//Not necessary for now
    //Setters
    void SetAlphaBasal(double alphaBasalInput){alphaBasal=alphaBasalInput;}
    //void SetAlphaStimmulusRest(double alphaStimmulusRest){alphaStimmulus=alphaStimmulusRest;}
    void SetAlphaExpDecay(double alphaStimExpDecay){alphaStimmulusExpDecay=alphaStimExpDecay;}
    // void SetNBasal(int nBasalInput){nBasal=nBasalInput;}
    // void SetKBasal(int kBasalInput){kBasal=kBasalInput;}
    // void SetKExponentialDecay(int kStimmulusExpDecayCalc){kStimmulusExpDecay=kStimmulusExpDecayCalc;}//Input is supposed to be exp(-dt/tau)
    // void SetNExponentialDecay(int nStimmulusExpDecayCalc){nStimmulusExpDecay=nStimmulusExpDecayCalc;}
    void SetMaxCount(int maxCountnew) {maxCount=maxCountnew;} //UNRESOLVED  Call in AllocateNewSynapse
    void SetPotentiationRatio(double ratio) {PotentiationDepressionRatio=ratio;}
    // void SetDepressionFlag(bool booleanFlag) {depressionFlagSTDP=booleanFlag;}
    // void SetPotentiationFlag(bool booleanFlag) {potentiationFlagSTDP=booleanFlag;}
    void SetUpdatedFlag(bool flag){updated=flag;}
    //Alpha methods
    void RecalculateAlphaResources();//Should be recalced at least once in the Spine Setup
    void DecayAlphaResources();
    void RecalcWeight(double weightResourceFactor);
    //Temp effets methods
    void AddTempResourcesToSpine(double alphaStimmulusInput);
    bool ApplyAllTemporaryEffectsOnPostspike(const DHashMap& STDPdecayMap);//input must be -1 if depression or STDPratio if potentiation, or the inverse swapping everything
    //void ApplyAllTempEffectsOnConflictPotentiation(double PotentiationDepressionRatio);//input must be -1 if depression or STDPratio if potentiation, or the inverse swapping everything
    bool ApplyAllTemporaryEffectsOnDepression(double expDecayFactor);//input must be -1 if depression or STDPratio if potentiation, or the inverse swapping everything

    //Stimmulus vector methods
    void TickStimmulusCounts();//Called in Reset(), but both should be mutually exclusive with AATE above
    void CullStimmulusVectors();//Called in Reset()
    //Profile methods
    std::valarray<double> GetIndividualSynapticProfile() const override;
    std::string GetIndividualSynapticProfileHeaderInfo() const override;
    //Virtual bools
    bool IgnoreCouplingStrength() override {return true;}
};

#endif