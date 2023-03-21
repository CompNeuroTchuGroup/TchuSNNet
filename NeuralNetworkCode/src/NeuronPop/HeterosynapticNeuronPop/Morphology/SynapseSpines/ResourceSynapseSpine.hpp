#ifndef _SYNAPSE_SPINE_RESOURCE_CLASS_HEADER_
#define _SYNAPSE_SPINE_RESOURCE_CLASS_HEADER_

#include "BranchedSynapseSpine.hpp"
#include <vector>

class ResourceSynapseSpine : public BranchedSynapseSpine {

    protected:
    //Resource model.
    double kBasal{1.0};
    double nBasal{1.0};

    double kStimmulus{1.0};
    double nStimmulus{1.0};

    double kExpDecay{1.0};
    double nExpDecay{1.0};

    double alphaResources{};

    std::vector<double> kStimmulusTempEffect{};//These need to be here for non-matrix management.
    std::vector<double> nStimmulusTempEffect{};
    std::vector<int> stimmulusEffectCount{};

    int maxCount{100};

    public:
    //Getters
    //int GetMaxCount(){return maxCount;}//Not necessary for now
    //Setters
    void SetMaxCount(int maxCountnew) {maxCount=maxCountnew;}
    //Alpha methods
    void RecalculateAlphaResources();//Should be recalced at least once in the Spine Setup
    void DecayAlphaResources();
    //Temp effets methods
    void AddTempResourcesToSpine(double kS, double nS);
    void ApplyAllTempEffects(bool potentiation);
    //Stimmulus vector methods
    void TickStimmulusCounts();
    void CullStimmulusVectors();
    //Profile methods
    std::valarray<double> GetIndividualSynapticProfile() const override;
    std::string GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif