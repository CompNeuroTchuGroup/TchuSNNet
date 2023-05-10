#ifndef _RESOURCE_HETERO_SYNAPTIC_PLASTICITY_STDP_BASED_HEADER_
#define _RESOURCE_HETERO_SYNAPTIC_PLASTICITY_STDP_BASED_HEADER_

#include "../BranchedMorphology.hpp"
#include "../SynapseSpines/ResourceSynapseSpine.hpp"
#include <numeric>
#include <unordered_map>
#include <map>
#include <unordered_set>

class BranchedMorphology;

typedef std::shared_ptr<ResourceBranch> RBranchPtr;
typedef std::shared_ptr<Branch> BranchPtr;
typedef std::shared_ptr<ResourceSynapseSpine> ResourceSpinePtr;
typedef std::shared_ptr<BranchedSynapseSpine> BranchedSpinePtr;
typedef std::shared_ptr<BaseSynapseSpine> BaseSpinePtr;
typedef std::unordered_map<int, double> DHashMap;
typedef std::unordered_map<int, DHashMap> NestedDHashMap;

class BranchedResourceHeteroSTDP : public BranchedMorphology {
//This class models a behaviour based on wi=beta*(alfai/(omega+sum(alfai))), where alfai represents the spine's resources as (Ks*expdt+Kbasal)/(Ns*expdt+Nbasal) with bumps on Ks and Ns
protected:
    //Synapse variables
    double alphaBasal{1.0};//LP and SP 
    double alphaStimmulusExpTau{1.0};//LP and SP store the tau
    double alphaStimmulusExpDecay{1.0}; //LP and SP has to already be calculated
    //Branch variables
    // size_t betaEventsWindowSize{};
    // int betaEventsPerTimestepThreshold{};
    double betaResourcePool{1.0}; //LP and SP
    // double betaUpTick{0.05};
    double omegaPassiveResourceOffset{1.0}; //LP and SP
    //Space-time kernel
    int kernelGapNumber{};//Obtained by dividing LP param by synGap and rounded to int
    int timeKernelLength{};//Obtained by dividing LP param by dt and rounded to int
    //The following two variables can be either ints or doubles. Doubles avoid rounding down and thus more precise area management
    double spaceTimeStepRelation{}; //Amount of timesteps per gap in kernel space. Allows iterating while keeping the kernel as a stair. 
    double kernelRadiusArea{};//Amount of timesteps from center of kernel to the border of kernel. 

    double timeKernelExpDecayConstant{1.0};//Strong decay needs small constants//LP and SP
    double spaceKernelExpDecayConstant{1.0};//LP and SP
    double ExpDecayConstantSTDP{1.0};//LP and SP

    NestedDHashMap kernelHashTable{};//This is more efficient than the map, but we need a hash function. TEST whether this is faster than vector<vector> or map<>
    DHashMap DecayHashTableSTDP{};//??Should STDP not decay? unordered_map? looks like vector is more ideal

    // double baseKStimmulusBump{1.0};//LP and SP
    // double baseNStimmulusBump{1.0};//LP and SP
    double baseAlphaStimmulusBump{1.0};//LP and SP

    //STDP and STDD parameters (MAKE SURE TO TUNE STDD, AS NEGATIVE WEIGHTS MUST BE AVOIDED)
    double PotentiationDepressionRatio{1.0};//Diagnose later negative weights //LP and SP (for now it can be ignored)

    //Max counters
    int MaxCountSTDP{};//This one is the lasting count for th spines //LP and SP?//dependent on dt?, default 10 ms assuming dt=1e-4
    //int branchMaxCountTrigger{};//LP and SP?//dependent on dt?, default 10 ms assuming dt=1e-4 //This variable should be equivalent to timeKernelLength

    //Counting
    int STDPDepressionCount{};//In relation to maxCountSTDP
    //Class object pointer vectors (virtual access)
    std::vector<ResourceSpinePtr> resourceSynapseData;
    std::vector<RBranchPtr> resourceBranches;

    //Record variables
    int totalLTPEvents{};
    int totalLTDEvents{};
    
public:

    //main Methods
    BranchedResourceHeteroSTDP()=default;
    explicit BranchedResourceHeteroSTDP(GlobalSimInfo * info);
    ~BranchedResourceHeteroSTDP() = default;
    
    void LoadParameters(std::vector<std::string> *input) override;//Remember to set all counts to maxCount    
    void SaveParameters(std::ofstream * stream, std::string neuronPreId) override;

    void SetUpBranchings(int remainingBranchingEvents, std::vector<int> anteriorBranches = std::vector<int>()) override;// Here we set up the vector with the branches
    void SetUpHashTables(); //Has to set up both time and space from the exp constants. Call in LP

    const std::string GetType() override {return str_BranchedResourceHeteroSTDP;};

    //Advect methods
    void advect() override;
    //Pairing functions
    void DetectPossiblePairing(RBranchPtr& branch);
    //bool CheckIfThereIsPairing(RBranchPtr& branch, int synapseIDinBranch);
    void SpaceTimeKernel(int branchSynapseID, int branchID);
    double CallKernelHashTable(int distanceToCenterInGaps, int timeDifference);
    //Plasticity events functions
    void ApplyEffects();//Here we increase the plasticity count of synapse and branch
    // void STDPPotentiation(ResourceSpinePtr& synapse);
    // void STDPDepression(ResourceSpinePtr& synapse);
    //Reset methods
    void Reset() override; //Wrapper plus clearing some of the vectors. Last Reset method to run in chronological order, where we call the ticks and the general upkeep
    void DeleteEffects();//Here, if counter==countMax, erase in that index the element of every vector (first store index, then REVERSE remove the removelist indexes with .rbegin and .rend instead of .begin and .end)
        //Container should be ordered by definition, but std::sort(array.begin(), array.end()) would ensure so.
    void TickAllCounts();//Last method called in Reset()
    void ClearSynapseSets();
    //Recalc methods. These methods have to be done per branch
    void RecalcAlphas(RBranchPtr& branch);//Run in LP
    void RecalcWeights(RBranchPtr& branch);//Run in LP
    void RecalcAlphaSums(RBranchPtr& branch);//Called inside recalc weights
    //Record methods
    void RecordPostSpike() override;
    void RecordExcitatoryPreSpike(int spikedSynapseId) override;//Here set the trigger count to 0

    //Allocation methods
    BaseSpinePtr AllocateNewSynapse(const HeteroCurrentSynapse& synapse) override; //Call the Branched one inside before setting all counters
        //Remember to set all counts to maxCount    
    //Record functions
    std::valarray<double> GetOverallSynapticProfile() override;
    std::string GetOverallSynapticProfileHeaderInfo() const override;
    //void CalcMorphoPlasticityEvents() override;
    //For debugging purposes
    void WeightDecay() override {throw;};//This function should never be called 
    void NormalizeWeights() override {return;};
};

#endif