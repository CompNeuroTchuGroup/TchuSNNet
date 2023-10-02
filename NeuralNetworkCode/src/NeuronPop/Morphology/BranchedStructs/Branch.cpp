#include "./Branch.hpp"
#include "Branch.hpp"

Branch::Branch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId):branchId{branchId},anteriorBranches{anteriorBranches},synapticGap{gap}, branchLength{branchLength}, branchSlots{static_cast<size_t>(std::round(branchLength/gap))}{//,morphoSynapseIDs(static_cast<size_t>(branchSlots), -1),branchSynapseIDs(static_cast<size_t>(branchSlots), -1), , spikedSyn(branchSlots, false) {
    //std::iota(uniqueSynapsePositionIDs.begin(),uniqueSynapsePositionIDs.end() , branchId*(branchSlots));
    //Constructor
}

Branch::Branch(double gap, double branchLength, int branchId):anteriorBranches{},branchId{branchId},synapticGap{gap}, branchLength{branchLength}, branchSlots{static_cast<size_t>(std::round(branchLength/gap))} {
}

void Branch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
    //Here we do all the function calls that could not be done in the constructor. 
    //This is done to adapt to the fact that synapses do not exist until ConnectNeurons() is called in the NeuralNetwork::Simulate() function.
    //For now, it is not necessary to do anything in this class
}
