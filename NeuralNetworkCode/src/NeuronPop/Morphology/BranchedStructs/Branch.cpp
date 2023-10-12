#include "./Branch.hpp"
#include "Branch.hpp"

Branch::Branch(std::vector<int> anteriorBranches,double gap, double branchLength,  int branchId):branchId{branchId},anteriorBranches{anteriorBranches},synapticGap{gap}, branchLength{branchLength}, branchSlots{static_cast<size_t>(std::round(branchLength/gap))}{//,morphoSynapseIDs(static_cast<size_t>(branchSlots), -1),branchSynapseIDs(static_cast<size_t>(branchSlots), -1), , spikedSyn(branchSlots, false) {
    //std::iota(uniqueSynapsePositionIDs.begin(),uniqueSynapsePositionIDs.end() , branchId*(branchSlots));
    //Constructor
}

Branch::Branch(double gap, double branchLength, int branchId):anteriorBranches{},branchId{branchId},synapticGap{gap}, branchLength{branchLength}, branchSlots{static_cast<size_t>(std::round(branchLength/gap))} {
}
