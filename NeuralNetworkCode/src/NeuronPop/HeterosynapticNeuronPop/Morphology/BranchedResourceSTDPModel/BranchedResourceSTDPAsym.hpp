#ifndef _SIMPLE_PLASTICITY_ONLY_BRANCHES_
#define _SIMPLE_PLASTICITY_ONLY_BRANCHES_

#include "../BranchedMorphology.hpp"
#include "../SynapseSpines/ResourceSynapseSpine.hpp"
#include <numeric>
#include <unordered_map>
#include <map>

class BranchedMorphology;

typedef std::unordered_map<int, double> DHashMap;
typedef std::unordered_map<int, DHashMap> SuperHashMap;

class BranchedResourceSTDPAsymmetric : public BranchedMorphology {
//This class models a behaviour based on wi=beta*(alfai/(omega+sum(alfai))), where alfai represents the spine's resources as (Ks*expdt+Kbasal)/(Ns*expdt+Nbasal) with bumps on Ks and Ns
protected:
    //Branch variables
    // size_t betaEventsWindowSize{};
    // int betaEventsPerTimestepThreshold{};
    double betaResourcePool{1.0};
    // double betaUpTick{0.05};
    double omegaPassiveResourceOffset{1.0};
    //Space-time kernel
    int kernelGapNumber{};//Obtained by dividing LP param by synGap and rounded to int
    int timeKernelLength{};//Obtained by dividing LP param by dt and rounded to int
    int spaceTimeStepRelation{};
    int kernelRadius{};

    double timeKernelDecayConstant{1.0};//Strong decay needs small constants//LP and SP
    double spaceKernelDecayConstant{1.0};//LP and SP
    double STDPDecayConstant{1.0};//LP and SP

    SuperHashMap kernelHashTable;//This is more efficient than the map, but we need a hash function. TEST whether this is faster than vector<vector> or map<>
    DHashMap STDPDecayHashTable;//??Should STDP not decay? unordered_map? looks like vector is more ideal

    // double baseKStimmulusBump{1.0};//LP and SP
    // double baseNStimmulusBump{1.0};//LP and SP
    double baseAlphaStimmulusBump{1.0};//LP and SP

    //STDP and STDD parameters (MAKE SURE TO TUNE STDD, AS NEGATIVE WEIGHTS MUST BE AVOIDED)
    double PotentiationDepressionRatio{1.00};//Diagnose later negative weights //LP and SP (for now it can be ignored)

    //Max counters
    int maxCount{100}; //dependent on dt?, default 10 ms assuming dt=1e-4 //LP and SP, this problem is centralized to this class now
    int& MaxCountSTDP = maxCount;//This one is the lasting count for th spines //LP and SP?
    int& branchMaxCountTrigger = maxCount;//LP and SP?

    //Counting
    int STDPDepressionCount{};//In relation to maxCountSTDP
    //Class object pointer vectors (virtual access)
    std::vector<std::shared_ptr<ResourceSynapseSpine>> synapseDataResources;
    std::vector<std::shared_ptr<ResourceBranch>> resourceBranches;
public:

    //main Methods
    BranchedResourceSTDPAsymmetric()=default;
    explicit BranchedResourceSTDPAsymmetric(GlobalSimInfo * info);
    ~BranchedResourceSTDPAsymmetric() = default;
    
    void LoadParameters(std::vector<std::string> *input) override;//Remember to set all counts to maxCount    
    void SaveParameters(std::ofstream * stream, std::string neuronPreId) override;

    void SetUpBranchings(int remainingBranchingEvents, std::vector<int> anteriorBranches = std::vector<int>()) override;// Here we set up the vector with the branches
    void SetUpHashTables(); //Has to set up both time and space from the exp constants. Call in LP

    const std::string GetType() override {return str_BranchedResourceSTDPAsymmetric;};

    
    //Advect methods
    void advect() override;
    void STDPPotentiation();
    void STDPDepression();
    void DetectPairing(std::vector<int> synapseIndexesToUpdate);
    void SetEffects(int synapseSpineIDinMorpho);//Here we calculate and store the effects of pairing (no matter the post state)
    void SpaceTimeKernel(int branchSynapseID, int branchID, int synapseSpineIDinMorpho);
    double CallKernelHashTable(int distanceToCenterInGaps, int timeDifference);
    void TickCounts();
    void ApplyEffects();//Here we increase the plasticity count of synapse and branch
    void DeleteEffects();//Here, if counter==countMax, erase in that index the element of every vector (first store index, then REVERSE remove the removelist indexes with .rbegin and .rend instead of .begin and .end)
        //Container should be ordered by definition, but std::sort(array.begin(), array.end()) would ensure so.
    void Reset() override; //Wrapper plus clearing some of the vectors. Last Reset method to run in chronological order, where we call the ticks and the general upkeep
    //These methods have to be done per branch
    void UpdateAlphaSum(std::shared_ptr<ResourceBranch> branch);//Run in LP
    void RecalcAlphas(std::shared_ptr<ResourceBranch> branch);//Run in LP
    void RecalcWeights(std::shared_ptr<ResourceBranch> branch);//Run in LP
    void RecalcAlphaSums(std::shared_ptr<ResourceBranch> branch);//Called inside recalc weights
    //Record methods
    void RecordPostSpike() override;
    void RecordExcitatoryPreSpike(unsigned long spikedSynapseId) override;//Here set the trigger count to 0

    //Getters
    std::valarray<double> GetIndividualSynapticProfile(unsigned long synapseId) const override;//Remember plasticity events per synapse and something else

    //Allocation methods
    std::shared_ptr<BaseSynapseSpine> AllocateNewSynapse(HeteroCurrentSynapse& synapse) override; //Call the Branched one inside before setting all counters
        //Remember to set all counts to maxCount    
    //For debugging purposes
    void WeightDecay() override {throw;};//This function should never be called 
    void NormalizeWeights() override {return;};
};

#endif