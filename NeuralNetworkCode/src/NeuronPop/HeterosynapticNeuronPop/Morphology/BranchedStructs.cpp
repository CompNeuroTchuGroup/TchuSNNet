#include "./BranchedStructs.hpp"
#include "BranchedStructs.hpp"

DendriticSubRegion::DendriticSubRegion(char regionID, std::vector<int> branchesInRegion): regionID{regionID}, branchesInRegion{branchesInRegion}
{
    //This is still under work
}

Branch::Branch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId):branchSlots{static_cast<size_t>(std::round(branchLength/gap))}, synapticGap{gap}, branchLength{branchLength}, anteriorBranches{anteriorBranches}, branchId{branchId}//,morphoSynapseIDs(static_cast<size_t>(branchSlots), -1),branchSynapseIDs(static_cast<size_t>(branchSlots), -1), , spikedSyn(branchSlots, false)
{
    //std::iota(uniqueSynapsePositionIDs.begin(),uniqueSynapsePositionIDs.end() , branchId*(branchSlots));
    //Constructor
}

void Branch::PostConnectSetUp(std::vector<BranchedSpinePtr> SynapseData)
{
    //Here we do all the function calls that could not be done in the constructor. 
    //This is done to adapt to the fact that synapses do not exist until ConnectNeurons() is called in the NeuralNetwork::Simulate() function.
    //For now, it is not necessary to do anything in this class
}

ResourceBranch::ResourceBranch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId, int triggerMaxCount): 
Branch(gap, branchLength, anteriorBranches, branchId), triggerCount(branchSlots, triggerMaxCount)//, potentiationCountSTDP(branchSlots, 10)//, maxCountSTDPPotentiation{branchMaxCountSTDPPotentiation}//, maxCountTrigger{branchMaxCountTrigger}//, plasticityEventsPerTimestepWindow(betaEventsWindowSize)
{
    //Constructor
    resouceBranchSynapseData= std::vector<ResourceSpinePtr>(branchSlots, nullptr);
}
void ResourceBranch::PostConnectSetUp(std::vector<BranchedSpinePtr> synapseData)
{
    //Here we do all the function calls that could not be done in the constructor. 
    //This is done to adapt to the fact that synapses do not exist until ConnectNeurons() is called in the NeuralNetwork::Simulate() function.
    SetUpSynapseData(synapseData);
}
void ResourceBranch::SetUpSynapseData(std::vector<BranchedSpinePtr> synapseData)
{
    //This function sets up the spine pointers in the branch, something not possible until ConnectNeurons() runs, thus its execution from PostConnectSetUp()
    for (BranchedSpinePtr synapse :synapseData){
        if (synapse->GetBranchId()==branchId){
            resouceBranchSynapseData.at(synapse->GetBranchPositionId())=std::dynamic_pointer_cast<ResourceSynapseSpine>(synapse);
        }
    }
    //std::sort(resouceBranchSynapseData.begin(), resouceBranchSynapseData.end(), BranchIDCompare); //This allows to do indexing of the synapse data using the branch ID
}

void ResourceBranch::TickAllCounts()
{
    //This function is supposed to increase the counts of all synapses
    for (int& count:this->triggerCount){
    //if (count<maxCount){
        count++;
    //}
    }
}

// void ResourceBranch::CheckIncreaseInBetaResources()
// {
//     //Assuming the plasticityEventsInThisTimestep is tracked properly:
//     plasticityEventsPerTimestepWindow.pop_front();
//     plasticityEventsPerTimestepWindow.push_back(plasticityBranchEventsThisTimestep);
//     if (plasticityBranchEventsThisTimestep>0){//This requirement is a logically coherent efficiency, as potentiation should not happen if there are no plasticity events.
//         plasticityBranchEventsTotal+=plasticityBranchEventsThisTimestep;
//         sumPlasticityEvents =std::accumulate(this->plasticityEventsPerTimestepWindow.begin(), this->plasticityEventsPerTimestepWindow.end(), 0);
//         if (sumPlasticityEvents>=plasticityEventsPerTimestepThreshold){
//             betaResourcePool += betaUpTick;
//             //here, do we reset the deque, or do we set a refractory period?
//         }
//         plasticityBranchEventsThisTimestep=0;
//     }
// }