#include "./AlphaBranch.hpp"

#include "AlphaBranch.hpp"

AlphaBranch::AlphaBranch(std::vector<int> anteriorBranches,double gap, double branchLength,  int branchId,
                         double preSynTraceDecay, double coopTraceDecay)
    : Branch(anteriorBranches,gap, branchLength,  branchId), spinePtrPosition(branchSlots, nullptr),
      preSynapticTraces(branchSlots), cooperativityTraces(branchSlots), preSynTraceDecay{preSynTraceDecay},
      coopTraceDecay{coopTraceDecay} {
}

AlphaBranch::AlphaBranch(double gap, double branchLength, int branchId, double preSynTraceDecay, double coopTraceDecay)
    : Branch(gap, branchLength, branchId), spinePtrPosition(branchSlots, nullptr), preSynapticTraces(branchSlots),
      cooperativityTraces(branchSlots), preSynTraceDecay{preSynTraceDecay}, coopTraceDecay{coopTraceDecay} {
}

void AlphaBranch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
    // Here we do all the function calls that could not be done in the
    // constructor. This is done to adapt to the fact that synapses do not
    // exist until ConnectNeurons() is called in the NeuralNetwork::Simulate()
    // function.
    for (BranchedSpinePtr synapse : spineData) {
        if (synapse->branchId == branchId) {
            spinePtrPosition.at(synapse->branchPositionId) = static_cast<AlphaSynapseSpine *>(synapse);
        }
    }
    // std::sort(resouceBranchSpineData.begin(), resouceBranchSpineData.end(),
    // BranchIDCompare); //This allows to do indexing of the synapse data using
    // the branch ID
}

void AlphaBranch::ApplyTracesOnSpinesLTP() {
    // This function exclusively performs LTP, as it is called during a
    // postspike. The cooperativity trace will be added to the presynaptic
    // trace to elicit heterosynaptic plasticity As it is trace-based, all
    // synapses are updated according to the trace. The idea behind the
    // equation is to make the trace heavily dependant on the presynaptic trace
    // as no glutamate means no calcium on postspike.
    for (size_t preSynIndex : std::ranges::views::iota(0u, branchSlots)) {
        if (spinePtrPosition.at(preSynIndex) != nullptr) {
            spinePtrPosition.at(preSynIndex)
                ->AddTraceToAlpha((1 + cooperativityTraces.at(preSynIndex)) *
                                  preSynapticTraces.at(preSynIndex)); // 1 can be substituted
                                                                      // with presynTrace
        }
    }
}

void AlphaBranch::ApplyTracesOnSpinesLTD(double postSynapticTrace, double biasLTD) {
    // This function exclusively performs LTD, as it is called during a
    // prespike. The cooperativity trace will be added to the presynaptic trace
    // to elicit heterosynaptic plasticity As it is trace-based, all synapses
    // are updated according to the trace. The idea behind the equation is to
    // make the trace heavily dependant on the presynaptic trace as no
    // glutamate means no calcium on postspike.
    for (int preSynIndex : spikedSpinesInTheBranch) {
        if (spinePtrPosition.at(preSynIndex) != nullptr) {
            spinePtrPosition.at(preSynIndex)
                ->AddTraceToAlpha(-(biasLTD - cooperativityTraces.at(preSynIndex)) *
                                  postSynapticTrace); // Negative to get LTD//1 can
                                                      // be substituted with
                                                      // postsynTrace
        }
    }
}

void AlphaBranch::DecayAllTraces() {
    // (std::placeholders::_1 is from std::placeholders, and represents future
    // arguments that will be passed to the std::function<T> object, which is
    // something equivalent to a lambda function).
    //  To understand what is going on :
    //  https://www.youtube.com/watch?v=1nNKeCFtBDI
    std::transform(std::execution::unseq, preSynapticTraces.begin(), preSynapticTraces.end(), preSynapticTraces.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, preSynTraceDecay));
    std::transform(std::execution::unseq, cooperativityTraces.begin(), cooperativityTraces.end(),
                   cooperativityTraces.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, coopTraceDecay));
}