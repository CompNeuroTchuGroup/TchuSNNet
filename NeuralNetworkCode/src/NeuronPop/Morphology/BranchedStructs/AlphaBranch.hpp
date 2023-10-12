//
// Created by Antoni Bertolin on 01.10.23
//
#ifndef _ALPHA_BRANCH_STRUCT_HEADER_
#define _ALPHA_BRANCH_STRUCT_HEADER_

#include "./Branch.hpp"
#include "../SynapseSpines/AlphaSynapseSpine.hpp"

struct AlphaBranch : public Branch {

    //Synapse access
    std::vector<AlphaSynapseSpine> alphaSpines;
    std::vector<AlphaSpinePtr> spinePtrPosition;//CAREFUL! THIS VECTOR IS NOT SORTED AT ANY POINT by SynapseBranchID
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
    AlphaBranch(std::vector<int>anteriorBranches,double gap, double branchLength,  int branchId, double preSynTraceDecay, double coopTraceDecay);
    AlphaBranch(double gap, double branchLength, int branchId, double preSynTraceDecay, double coopTraceDecay);    
    void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) override;
        //Count related functions
    void ApplyTracesOnSpinesLTP();
    void ApplyTracesOnSpinesLTD(double postSynapticTrace, double biasLTD);
    void DecayAllTraces();//Use ternary operator. Called in Reset()
    //void TickCounts(std::vector<int>& countVector);
    //void CheckIncreaseInBetaResources(); //Here I have to add current, delete last, sum and check against threshold. Called in Reset()
    //void UpdateAlphaSum(); //This dunction cannot be implemented in this struct as the struct has no access to the spine pointers
    //void IncreasePlasticityCounter() override {plasticityBranchEventsThisTimestep++;plasticityBranchEventsTotal++;}
};
#endif