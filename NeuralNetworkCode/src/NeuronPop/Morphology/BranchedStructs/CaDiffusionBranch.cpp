#include "CaDiffusionBranch.hpp"

CaDiffusionBranch::CaDiffusionBranch(std::vector<int> anteriorBranches, double gap, double branchLength, int branchId,
                                     cadouble initialCalcium, double initialResources, double initialWeight,
                                     double kinasesTotal, double phosphatasesTotal, double reaction1Ctt,
                                     double reaction2Ctt, double reaction3Ctt, double reaction4Ctt, double reaction5Ctt,
                                     double reaction6Ctt, cadouble caDiffusionFct, double resourceDiffusionFct,
                                     cadouble caDecayFct, int prespikeDelaySteps,double camSaturationSteep, double camSaturationInflection)
    : Branch(anteriorBranches, gap, branchLength, branchId), kinasesTotal{kinasesTotal},
      phosphatasesTotal{phosphatasesTotal}, reaction1Ctt{reaction1Ctt}, reaction2Ctt{reaction2Ctt},
      reaction3Ctt{reaction3Ctt}, reaction4Ctt{reaction4Ctt}, reaction5Ctt{reaction5Ctt}, reaction6Ctt{reaction6Ctt},
      caDiffusionFct{caDecayFct / std::pow(gap, 2)}, resourceDiffusionFct{resourceDiffusionFct / std::pow(gap, 2)},
      caDecayFct{caDecayFct}, camSaturationSteep{camSaturationSteep}, camSaturationInflection{camSaturationInflection} {

    CaResSpines.resize(branchLength / gap, CaResSynapseSpine(kinasesTotal, phosphatasesTotal));
    for (CaResSynapseSpine& spine : CaResSpines){
      spine.resourcesAvailable=initialResources;
      spine.calciumFree=initialCalcium;
      spine.weight=initialWeight;
    }
    waitingMatrix.resize(prespikeDelaySteps, std::vector<cadouble>(branchLength / gap, 0));
}

CaDiffusionBranch::CaDiffusionBranch(double gap, double branchLength, int branchId, cadouble initialCalcium,
                                     double initialResources, double initialWeight, double kinasesTotal,
                                     double phosphatasesTotal, double reaction1Ctt, double reaction2Ctt,
                                     double reaction3Ctt, double reaction4Ctt, double reaction5Ctt, double reaction6Ctt,
                                     cadouble caDiffusionFct, double resourceDiffusionFct, cadouble caDecayFct,
                                     int prespikeDelaySteps,double camSaturationSteep, double camSaturationInflection)
    : Branch(gap, branchLength, branchId), kinasesTotal{kinasesTotal}, phosphatasesTotal{phosphatasesTotal},
      reaction1Ctt{reaction1Ctt}, reaction2Ctt{reaction2Ctt}, reaction3Ctt{reaction3Ctt}, reaction4Ctt{reaction4Ctt},
      reaction5Ctt{reaction5Ctt}, reaction6Ctt{reaction6Ctt}, caDiffusionFct{caDecayFct / std::pow(gap, 2)},
      resourceDiffusionFct{resourceDiffusionFct / std::pow(gap, 2)}, caDecayFct{caDecayFct}, camSaturationSteep{camSaturationSteep}, camSaturationInflection{camSaturationInflection} {

    CaResSpines.resize(branchLength / gap, CaResSynapseSpine(kinasesTotal, phosphatasesTotal));
    for (CaResSynapseSpine& spine : CaResSpines){
      spine.resourcesAvailable=initialResources;
      spine.calciumFree=initialCalcium;
      spine.weight=initialWeight;
    }
    waitingMatrix.resize(prespikeDelaySteps, std::vector<cadouble>(branchLength / gap, 0));
}

void CaDiffusionBranch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
}

void CaDiffusionBranch::PreSpikeCalciumInflux(TStepInt step) {
    for (size_t index : std::ranges::views::iota(0u, CaResSpines.size())) {
        CaResSpines.at(index).calciumFree += waitingMatrix.at(step).at(index);
    }
    // for (size_t index : std::ranges::views::iota(0u, prespikeWaitingMatrix.size()-1)){
    //   prespikeWaitingMatrix.at(index).swap(prespikeWaitingMatrix.at(index+1));
    // }
    std::fill(PAR_UNSEQ,waitingMatrix.at(step).begin(), waitingMatrix.at(step).end(), 0.0);
}

// void CaDiffusionBranch::PreSpikeCalciumInflux() {
//   for (size_t index : std::ranges::views::iota(0u, CaDiffSpines.size())){
//     CaDiffSpines.at(index).calciumFree += prespikeWaitingMatrix.at(0).at(index);
//   }
//   for (size_t index : std::ranges::views::iota(0u, prespikeWaitingMatrix.size()-1)){
//     prespikeWaitingMatrix.at(index).swap(prespikeWaitingMatrix.at(index+1));
//   }
//   std::fill(prespikeWaitingMatrix.back().begin(), prespikeWaitingMatrix.back().end(),0.0);
// }

void CaDiffusionBranch::ExecuteReactions() {
    // I apologize in advance, as this function requires being run with spatial and temporal locality for max speed.
    // Spatial locality is achieved in the spine vector (not pointers). Temporal locality is achieved by making a mess
    // of a function
    std::for_each(CaResSpines.begin(), CaResSpines.end(), [](CaResSynapseSpine &spine) {
        spine.RestoreSpine();
    }); // NOT IN A FOR EACH. ALL REACTIONS ARE EXECUTED SEQUENTIALLY
    // Both diffusion functions are executed simultaneously for temporal locality
}

void CaDiffusionBranch::CalciumCheck() {
    if (std::any_of(PAR_UNSEQ,CaResSpines.begin(), CaResSpines.end(),
                    [](CaResSynapseSpine &spine) { return spine.calciumFree < 0; })) {
        throw "Calcium got negative";
    }
}