#include "./AlphaBranch.hpp"

#include "AlphaBranch.hpp"

AlphaBranch::AlphaBranch(std::vector<int> anteriorBranches, double gap, double branchLength, int branchId, double preSynTraceDecay,
                         double coopTraceDecay, double alphaExpDecay, double alphaStimBump, double alphaBasal, double betaResourcePool, double omegaOffset)
    : Branch(anteriorBranches, gap, branchLength, branchId), spinePtrPosition(branchSlots, nullptr), preSynapticTraces(branchSlots),
      cooperativityTraces(branchSlots), preSynTraceDecay{preSynTraceDecay}, coopTraceDecay{coopTraceDecay}, alphaExpDecay{alphaExpDecay},
      alphaStimBump{alphaStimBump}, alphaBasal{alphaBasal}, betaResourcePool{betaResourcePool},omegaOffset{omegaOffset} {
}

AlphaBranch::AlphaBranch(double gap, double branchLength, int branchId, double preSynTraceDecay, double coopTraceDecay, double alphaExpDecay,
                         double alphaStimBump, double alphaBasal, double betaResourcePool, double omegaOffset)
    : Branch(gap, branchLength, branchId), spinePtrPosition(branchSlots, nullptr), preSynapticTraces(branchSlots), cooperativityTraces(branchSlots),
      spineAlphaStims(branchSlots), tracePlaceholder(branchSlots, 1.0), preSynTraceDecay{preSynTraceDecay}, coopTraceDecay{coopTraceDecay},
      alphaExpDecay{alphaExpDecay}, alphaStimBump{alphaStimBump}, alphaBasal{alphaBasal}, betaResourcePool{betaResourcePool},omegaOffset{omegaOffset} {
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
void AlphaBranch::DecayAlphaResources() {
    std::transform(PAR_UNSEQ, spineAlphaStims.begin(), spineAlphaStims.end(), spineAlphaStims.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, alphaExpDecay));
}

void AlphaBranch::ComputeAlphaResources() {
    for (size_t preSynIndex : std::ranges::views::iota0ull, branchSlots)) {
        if (spinePtrPosition.at(preSynIndex) != nullptr) {
            spinePtrPosition.at(preSynIndex)->alphaResources = alphaBasal + spineAlphaStims.at(preSynIndex);
            // We append the decay here, its faster? Or not
        }
    }
    DecayAlphaResources();
}

void AlphaBranch::ComputeAlphaSums() {
    ComputeAlphaResources();
    alphaTotalSum = std::accumulate(alphaSpines.begin(), alphaSpines.end(), 0.0,
                                    [](double accumulator, const AlphaSynapseSpine &const spine) { return accumulator + spine.alphaResources; });
}

void AlphaBranch::ComputeWeights() {
    ComputeAlphaSums();
    double resourceFactor{betaResourcePool /(omegaOffset + alphaTotalSum)};
    std::for_each(PAR_UNSEQ, alphaSpines.begin(), alphaSpines.end(),[resourceFactor](AlphaSynapseSpine& spine){
        spine.weight=spine.alphaResources*resourceFactor;
    });
}

void AlphaBranch::ApplyTracesOnSpinesLTP() {
    spikedSpinesInTheBranch.clear();
    // This function exclusively performs LTP, as it is called during a
    // postspike. The cooperativity trace will be added to the presynaptic
    // trace to elicit heterosynaptic plasticity As it is trace-based, all
    // synapses are updated according to the trace. The idea behind the
    // equation is to make the trace heavily dependant on the presynaptic trace
    // as no glutamate means no calcium on postspike.
    // for (size_t preSynIndex : std::ranges::views::iota0ull, branchSlots)) {
    //     if (spinePtrPosition.at(preSynIndex) != nullptr) {
    //         spinePtrPosition.at(preSynIndex)
    //             ->AddTraceToAlpha(((1 + cooperativityTraces.at(preSynIndex)) *
    //                               preSynapticTraces.at(preSynIndex))); // 1 can be substituted
    //                                                                   // with presynTrace

    //     }
    // }
    std::fill(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), 1.0);
    std::transform(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), cooperativityTraces.begin(), tracePlaceholder.begin(),
                   std::plus<double>());
    std::transform(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), preSynapticTraces.begin(), tracePlaceholder.begin(),
                   std::multiplies<double>());
    std::transform(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), tracePlaceholder.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, alphaStimBump));
    std::transform(PAR_UNSEQ, spineAlphaStims.begin(), spineAlphaStims.end(), tracePlaceholder.begin(), spineAlphaStims.begin(), std::plus<double>());
}

void AlphaBranch::ApplyTracesOnSpinesLTD(double postSynapticTrace, double biasLTD) {
    // This function exclusively performs LTD, as it is called during a
    // prespike. The cooperativity trace will be added to the presynaptic trace
    // to elicit heterosynaptic plasticity As it is trace-based, all synapses
    // are updated according to the trace. The idea behind the equation is to
    // make the trace heavily dependant on the presynaptic trace as no
    // glutamate means no calcium on postspike.
    // for (int preSynIndex : spikedSpinesInTheBranch) {
    //     if (spinePtrPosition.at(preSynIndex) != nullptr) {
    //         spinePtrPosition.at(preSynIndex)
    //             ->AddTraceToAlpha(-(biasLTD - cooperativityTraces.at(preSynIndex)) *
    //                               postSynapticTrace); // Negative to get LTD//1 can
    //                                                   // be substituted with
    //                                                   // postsynTrace
    //     }
    // }
    std::fill(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), biasLTD);
    std::transform(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), cooperativityTraces.begin(), tracePlaceholder.begin(),
                   std::minus<double>());
    std::transform(PAR_UNSEQ, tracePlaceholder.begin(), tracePlaceholder.end(), tracePlaceholder.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, alphaStimBump * postSynapticTrace));
    std::transform(PAR_UNSEQ, spineAlphaStims.begin(), spineAlphaStims.end(), tracePlaceholder.begin(), spineAlphaStims.begin(),
                   [this](double alphaStim, double traceAdd) { return std::max(alphaStim - traceAdd, -this->alphaBasal); });
}

void AlphaBranch::DecayAllTraces() {
    // (std::placeholders::_1 is from std::placeholders, and represents future
    // arguments that will be passed to the std::function<T> object, which is
    // something equivalent to a lambda function).
    //  To understand what is going on :
    //  https://www.youtube.com/watch?v=1nNKeCFtBDI
    std::transform(std::execution::unseq, preSynapticTraces.begin(), preSynapticTraces.end(), preSynapticTraces.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, preSynTraceDecay));
    std::transform(std::execution::unseq, cooperativityTraces.begin(), cooperativityTraces.end(), cooperativityTraces.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, coopTraceDecay));
}