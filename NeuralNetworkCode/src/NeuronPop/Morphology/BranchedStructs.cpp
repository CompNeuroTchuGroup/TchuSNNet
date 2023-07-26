#include "./BranchedStructs.hpp"
#include "BranchedStructs.hpp"

Branch::Branch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId):branchId{branchId},anteriorBranches{anteriorBranches},synapticGap{gap}, branchLength{branchLength}, branchSlots{static_cast<size_t>(std::round(branchLength/gap))}{//,morphoSynapseIDs(static_cast<size_t>(branchSlots), -1),branchSynapseIDs(static_cast<size_t>(branchSlots), -1), , spikedSyn(branchSlots, false) {
    //std::iota(uniqueSynapsePositionIDs.begin(),uniqueSynapsePositionIDs.end() , branchId*(branchSlots));
    //Constructor
    spikedSpinesInTheBranch.reserve(branchSlots);
}

void Branch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
    //Here we do all the function calls that could not be done in the constructor. 
    //This is done to adapt to the fact that synapses do not exist until ConnectNeurons() is called in the NeuralNetwork::Simulate() function.
    //For now, it is not necessary to do anything in this class
}

ResourceTraceBranch::ResourceTraceBranch(double gap, double branchLength, std::vector<int> anteriorBranches, int branchId, double preSynTraceDecay, double coopTraceDecay): Branch(gap, branchLength, anteriorBranches, branchId),rBranchSpineData(branchSlots, nullptr),
preSynapticTraces(branchSlots),cooperativityTraces(branchSlots), preSynTraceDecay{preSynTraceDecay}, coopTraceDecay{coopTraceDecay}{
}
void ResourceTraceBranch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
    //Here we do all the function calls that could not be done in the constructor. 
    //This is done to adapt to the fact that synapses do not exist until ConnectNeurons() is called in the NeuralNetwork::Simulate() function.
    for (BranchedSpinePtr synapse : spineData){
        if (synapse->GetBranchId()==branchId){
            rBranchSpineData.at(synapse->GetBranchPositionId())=static_cast<ResourceSynapseSpine*>(synapse);
        }
    }
    //std::sort(resouceBranchSpineData.begin(), resouceBranchSpineData.end(), BranchIDCompare); //This allows to do indexing of the synapse data using the branch ID
}

void ResourceTraceBranch::ApplyTracesOnSpinesLTP(double potDepRatio) {
    //This function exclusively performs LTP, as it is called during a postspike. The cooperativity trace will be added to the presynaptic trace to elicit heterosynaptic plasticity
    //As it is trace-based, all synapses are updated according to the trace.
    //The idea behind the equation is to make the trace heavily dependant on the presynaptic trace as no glutamate means no calcium on postspike.
    for (size_t preSynIndex : std::ranges::views::iota(0u, branchSlots)){
        if (rBranchSpineData.at(preSynIndex)!=nullptr){
            rBranchSpineData.at(preSynIndex)->AddTraceToAlpha((1+cooperativityTraces.at(preSynIndex))*preSynapticTraces.at(preSynIndex)*potDepRatio);//1 can be substituted with presynTrace
        }
    }
}

void ResourceTraceBranch::ApplyTracesOnSpinesLTD(double postSynapticTrace) {
    //This function exclusively performs LTD, as it is called during a prespike. The cooperativity trace will be added to the presynaptic trace to elicit heterosynaptic plasticity
    //As it is trace-based, all synapses are updated according to the trace.
    //The idea behind the equation is to make the trace heavily dependant on the presynaptic trace as no glutamate means no calcium on postspike.
    for (int preSynIndex : spikedSpinesInTheBranch){
        if (rBranchSpineData.at(preSynIndex)!=nullptr){
            rBranchSpineData.at(preSynIndex)->AddTraceToAlpha(-(1-cooperativityTraces.at(preSynIndex))*postSynapticTrace);//Negative to get LTD//1 can be substituted with postsynTrace
        }
    }
}

void ResourceTraceBranch::DecayAllTraces() {
    // (std::placeholders::_1 is from std::placeholders, and represents future arguments that will be passed to the std::function<T> object, which is something equivalent to a lambda function).
    //  To understand what is going on : https://www.youtube.com/watch?v=1nNKeCFtBDI
    std::transform(std::execution::unseq,preSynapticTraces.begin(), preSynapticTraces.end(), preSynapticTraces.begin(),std::bind(std::multiplies<double>(),std::placeholders::_1, preSynTraceDecay));
    std::transform(std::execution::unseq,cooperativityTraces.begin(), cooperativityTraces.end(), cooperativityTraces.begin(),std::bind(std::multiplies<double>(),std::placeholders::_1, coopTraceDecay));
}