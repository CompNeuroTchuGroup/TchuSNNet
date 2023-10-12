//
// Created by Antoni Bertolin on 15.09.23
//
#ifndef RESOURCE_HETERO_SYNAPTIC_PLASTICITY_SHARED_RESOURCES_STDP_TRACE_BASED_HEADER_
#define RESOURCE_HETERO_SYNAPTIC_PLASTICITY_SHARED_RESOURCES_STDP_TRACE_BASED_HEADER_

//List of forward declarations needed to break circular dependencies
// class BranchedMorphology;
struct BranchTargeting;
#include "../BranchedStructs/AlphaBranch.hpp"
#include "../BranchedMorphology.hpp"
#include "../../../GlobalFunctions.hpp"
#include <numeric>

// DEPRECATED
class AlphaRSharedBHSTDP : public BranchedMorphology {
//This class models a behaviour based on wi=beta*(alfai/(omega+sum(alfai))), where alfai represents the spine's resources as (Ks*expdt+Kbasal)/(Ns*expdt+Nbasal) with bumps on Ks and Ns
//It is mathematically provable that this model is equivlaent to the original in a single branch.
protected:
    //Synapse variables
    //Synapse variables
    double alphaBasal{1.0};//LP and SP 
    double alphaStimulusTau{1.0};//LP and SP store the tau
    double alphaStimulusExpDecay{1.0}; //LP and SP has to already be calculated
    //Branch variables
    // size_t betaEventsWindowSize{};
    // int betaEventsPerTimestepThreshold{};
    double betaResourcePool{1.0}; //LP and SP
    double alphaSum{1.0};
    double resourceMultiplier{1.0};
    // double betaUpTick{0.05};
    double omegaOffset{1.0}; //LP and SP

    double cooperativityTau{1.0};//Strong decay needs small constants//LP and SP
    double coopExpDecay{1.0};//Strong decay needs small constants//LP and SP
    double spaceProfileLambda{1.0};//LP and SP
    double tauSTDP{1.0};//LP and SP
    double STDPExpDecay{1.0};

    std::vector<double> spatialProfile{};//This is more efficient than the map, but we need a hash function. TEST whether this is faster than vector<vector> or map<>
    // std::unordered_map<int, double> decayHashTableSTDP{};//??Should STDP not decay? unordered_map? looks like vector is more ideal

    // double baseKStimulusBump{1.0};//LP and SP
    // double baseNStimulusBump{1.0};//LP and SP
    double baseAlphaStimBump{1.0};//LP and SP

    //STDP and STDD parameters (MAKE SURE TO TUNE STDD, AS NEGATIVE WEIGHTS MUST BE AVOIDED)
    double biasLTD{1.0};

    //PostSynaptic trace
    double postSynapticTrace{};
    std::vector<AlphaSpinePtr> resourceSpineData;
    std::vector<AlphaBranchPtr> alphaBranches;
    
public:

    //main Methods
    AlphaRSharedBHSTDP()=default;
    explicit AlphaRSharedBHSTDP(GlobalSimInfo* infoGlobal);
    ~AlphaRSharedBHSTDP() override = default;
    
    void LoadParameters(const std::vector<FileEntry>& morphologyParameters) override;
    void CheckParameters(const std::vector<FileEntry>& parameters) override;    
    void SaveParameters(std::ofstream& wParameterStream, std::string neuronIdentificator) const override;

    int CreateBranch(std::vector<int> anteriorBranches) override;
    void SetUpHashTable(); //Has to set up both time and space from the exp constants. Call in LP
    
    std::string GetType() const override {return IDstringTraceResourceHSTDP;};

    //Advect methods
    void Advect() override;
    //Pairing functions
    // void UpdateCoopTrace(const ResourceTraceBranch* const branch);
    bool CheckIfPreSpikeHappened();
    //bool CheckIfThereIsPairing(RBranchPtr branch, int synapseIDinBranch);
    void ApplyCoopTraceSpatialProfile(int branchSpineID, AlphaBranch* const branchID);
    // double CallKernelHashTable(int distanceToCenterInGaps);
    //Plasticity events functions
    // void ApplyEffects();//Here we increase the plasticity count of synapse and branch
    double GetPostSynapticTrace(){return postSynapticTrace;}
    double GetLTDBias(){return biasLTD;}
    // void STDPPotentiation(ResourceSpinePtr synapse);
    // void STDPDepression(ResourceSpinePtr synapse);
    //Reset methods
    void Reset() override; //Wrapper plus clearing some of the vectors. Last Reset method to run in chronological order, where we call the ticks and the general upkeep
    // void DeleteEffects();//Here, if counter==countMax, erase in that index the element of every vector (first store index, then REVERSE remove the removelist indexes with .rbegin and .rend instead of .begin and .end)
        //Container should be ordered by definition, but std::sort(array.begin(), array.end()) would ensure so.
    void DecayAllTraces();//Last method called in Reset()
    // void ClearSynapseSets();
    //Recalc methods. These methods have to be done per branch
    void ComputeAlphas();//Run in LP
    void ComputeWeights();//Run in LP
    void ComputeAlphaSums();//Called inside recalc weights
    //Record methods
    void RecordPostSpike() override;
    void RecordExcitatoryPreSpike(int spikedSpineId) override;//Here set the trigger count to 0
    void PostConnectSetUp() override;
    //Allocation methods
    BaseSpinePtr AllocateNewSynapse(const BranchTargeting& synapse) override; //Call the Branched one inside before setting all counters
        //Remember to set all counts to maxCount    
    //Record functions
    std::vector<double> GetOverallSynapticProfile() const override;
    std::string GetOverallSynapticProfileHeaderInfo() const override;
    //void CalcMorphoPlasticityEvents() override;
    //For debugging purposes
    bool IgnoreJDParameters() const override {return true;}
    void WeightDecay() override {throw "Unintended call of AlphaResourceHSTDP::WeightDecay";};//This function should never be called 
    void NormalizeWeights() override {return;};
};

#endif