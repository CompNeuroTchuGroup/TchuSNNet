#ifndef _BRANCHED_STRUCTS_HEADER_
#define _BRANCHED_STRUCTS_HEADER_

#include "./SynapseSpines/ResourceSynapseSpine.hpp"
#include <vector>
#include <deque>
#include <unordered_set>
#include <numeric>
#include <memory>
#include <algorithm>

typedef std::unordered_map<int, int> IHashMap;
typedef std::shared_ptr<ResourceSynapseSpine> ResourceSpinePtr;
typedef std::shared_ptr<BranchedSynapseSpine> BranchedSpinePtr;

struct DendriticSubRegion{ //Still under work
    const char regionID;
    const std::vector<int> branchesInRegion; //I will have to read this from the morphology LP, every DendriticSubRegion is a line, first input is ID, rest is branchIDs. Then in Synapse you put the DendriticSubRegion where the synapse goes. 
    DendriticSubRegion(char regionID, std::vector<int> branchesInRegion);
};

struct BranchTargeting{ //This is esentially a wrapper for HCS different targeting strategies of postsynaptic branches
    int targetBranch{};
    bool setTargetBranch{false};
    bool randomTargetBranch{false};
    bool orderedTargetBranch{false};
    char DendriticSubRegion{'0'};
};

struct Branch{
    //ID
    const int branchId{};
    //Branched vars
    const std::vector<int> anteriorBranches{}; 

    const double synapticGap{};    //For now these are identical to the morphology ones, but we will see in the future
    const double branchLength{};

    size_t branchSlots{};

    std::deque<int> openSynapsesSlots{};
    //For the actual checks
    // std::vector<bool> spikedSyn{};
    std::vector<int> spikedSynapsesInTheBranch{};
    std::vector<int> synapseSlotClosedIndex{};
    //std::vector<int> morphoSynapseIDs{};//This data variable is no longer relevant
    IHashMap synapseSlotToMorphoIDMap{};
    //Resource paradigm variables
    //int plasticityBranchEventsTotal{};
    int LTDevents{};
    int LTPevents{};
    //Methods
    //Branch()=default;
    Branch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId);
    virtual void PostConnectSetUp(std::vector<BranchedSpinePtr> synapseData);
    //virtual void IncreasePlasticityCounter(){plasticityBranchEventsTotal++;}
    void IncreasePotentiationCount(){LTPevents++;}
    void IncreaseDepressionCount(){LTDevents++;}
};

struct ResourceBranch : public Branch {

    //Synapse access
    std::vector<ResourceSpinePtr> resouceBranchSynapseData;//CAREFUL! THIS VECTOR IS NOT SORTED AT ANY POINT by SynapseBranchID

    //std::set<int> updatedSynapseSpines{};//Removed due to the fact that pairing in the same timestep should have a greater effect than when paired with greater tiem difference
    std::vector<int> updatedAlphaEffects{};//IDs are branch slots, Used in depression, cleared every timestep
    // int plasticityEventsPerTimestepThreshold{};
    //int sumPlasticityEvents{};//Theoretically temporary

    //double betaResourcePool{1.0}; //We use the one in the Morpho object
    // double betaUpTick{0.05};
    // /double omegaPassiveResourceOffset{1.0}; //We use the one in the Morpho object
    double alphaTotalSum{};
    double weightResourceFactor{};
    //Counts (all of these are not general, which is kind of an issue for generality)
    //int maxCountSTDPPotentiation{};
    //int maxCountTrigger{};
    std::vector<int> triggerCount{};//Here is where we look for counts under 10
    //std::vector<int> potentiationCountSTDP{};//Size of these has to be the amount of synapses in the branch (equal to triggerCount, length divided by gaps) //REMOVED due to redundancy
    //Misc
    // int plasticityBranchEventsThisTimestep{};//Has to be set to zero
    // std::deque<int> plasticityEventsPerTimestepWindow{};//Every timestep I push_front the events in the timestep and delete the last element with pop or erase. Then sum and check if over the threshold.
        //Do we clear() after a trigger or not? Most of the function would be the same, or put some refractory period UNRESOLVED
        //Here we could create a false history of plasticity events
    //Methods
    ResourceBranch(double gap, double branchLength, std::vector<int>anteriorBranches, int branchId, int triggerMaxCount);
    void PostConnectSetUp(std::vector<BranchedSpinePtr> synapseData) override;
    void SetUpSynapseData(std::vector<BranchedSpinePtr> synapseData);
        //Count related functions
    void TickAllCounts();//Use ternary operator. Called in Reset()
    //void TickCounts(std::vector<int>& countVector);
    //void CheckIncreaseInBetaResources(); //Here I have to add current, delete last, sum and check against threshold. Called in Reset()
    //void UpdateAlphaSum(); //This dunction cannot be implemented in this struct as the struct has no access to the spine pointers
    //void IncreasePlasticityCounter() override {plasticityBranchEventsThisTimestep++;plasticityBranchEventsTotal++;}
};
#endif