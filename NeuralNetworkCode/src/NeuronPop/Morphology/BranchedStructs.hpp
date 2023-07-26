//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _BRANCHED_STRUCTS_HEADER_
#define _BRANCHED_STRUCTS_HEADER_

#include "./SynapseSpines/ResourceSynapseSpine.hpp"
#include <vector>
#include <deque>
#include <unordered_set>
#include <numeric>
#include <memory>
#include <algorithm>
#include <functional>
#include <execution>

struct Branch{
    //ID
    const int branchId{};
    //Branched vars
    const std::vector<int> anteriorBranches{}; 

    const double synapticGap{};    //For now these are identical to the morphology ones, but we will see in the future
    const double branchLength{};

    size_t branchSlots{};

    std::deque<int> openSpineSlots{};
    //For the actual checks
    // std::vector<bool> spikedSyn{};
    std::vector<int> spikedSpinesInTheBranch{};
    std::vector<int> synapseSlotClosedIndex{};
    //std::vector<int> morphoSynapseIDs{};//This data variable is no longer relevant
    //IHashMap synapseSlotToMorphoIDMap{};
    //Resource paradigm variables
    //int plasticityBranchEventsTotal{};
    //Methods
    //Branch()=default;
    Branch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId);
    virtual void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData);
    //virtual void IncreasePlasticityCounter(){plasticityBranchEventsTotal++;}
};

struct ResourceTraceBranch : public Branch {

    //Synapse access
    std::vector<ResourceSynapseSpine*> rBranchSpineData;//CAREFUL! THIS VECTOR IS NOT SORTED AT ANY POINT by SynapseBranchID
    double alphaTotalSum{};
    double resourceFactor{};
    std::vector<double> preSynapticTraces;//Here is where we look for counts under 10
    std::vector<double> cooperativityTraces;//Here is where we look for counts under 10
    double preSynTraceDecay{};
    double coopTraceDecay{};
    //std::vector<int> potentiationCountSTDP{};//Size of these has to be the amount of synapses in the branch (equal to triggerCount, length divided by gaps) //REMOVED due to redundancy
    //Misc
    // int plasticityBranchEventsThisTimestep{};//Has to be set to zero
    // std::deque<int> plasticityEventsPerTimestepWindow{};//Every timestep I push_front the events in the timestep and delete the last element with neuronPop or erase. Then sum and check if over the threshold.
        //Do we clear() after a trigger or not? Most of the function would be the same, or put some refractory period UNRESOLVED
        //Here we could create a false history of plasticity events
    //Methods
    ResourceTraceBranch(double gap, double branchLength, std::vector<int>anteriorBranches, int branchId, double preSynTraceDecay, double coopTraceDecay);
    void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) override;
        //Count related functions
    void ApplyTracesOnSpinesLTP(double potDepRatio);
    void ApplyTracesOnSpinesLTD(double postSynapticTrace);
    void DecayAllTraces();//Use ternary operator. Called in Reset()
    //void TickCounts(std::vector<int>& countVector);
    //void CheckIncreaseInBetaResources(); //Here I have to add current, delete last, sum and check against threshold. Called in Reset()
    //void UpdateAlphaSum(); //This dunction cannot be implemented in this struct as the struct has no access to the spine pointers
    //void IncreasePlasticityCounter() override {plasticityBranchEventsThisTimestep++;plasticityBranchEventsTotal++;}
};
#endif