//
// Created by Antoni Bertolin on 01.10.23
//
#ifndef _BRANCH_STRUCT_HEADER_
#define _BRANCH_STRUCT_HEADER_

#include "../SynapseSpines/ResourceSynapseSpine.hpp"
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
    Branch(double gap, double branchLength, int branchId);
    virtual void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData);
    //virtual void IncreasePlasticityCounter(){plasticityBranchEventsTotal++;}
};

#endif