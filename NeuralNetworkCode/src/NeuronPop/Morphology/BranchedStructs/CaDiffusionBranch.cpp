#include "CaDiffusionBranch.hpp"

CaDiffusionBranch::CaDiffusionBranch(std::vector<int> anteriorBranches, double gap, double branchLength, int branchId, Constants constants)
    : Branch(anteriorBranches, gap, branchLength, branchId), constants{constants} {
    CaResSpines.resize(branchLength / gap, CaResSynapseSpine());
    for (CaResSynapseSpine &spine : CaResSpines) {
        spine.resourcesAvailable = constants.initialResources;
        spine.weight             = constants.initialWeight;
    }
    // waitingMatrix.resize(prespikeDelaySteps, std::vector<double>(branchLength / gap, constants.calciumInfluxBasal));
}

CaDiffusionBranch::CaDiffusionBranch(double gap, double branchLength, int branchId, Constants constants) : Branch(gap, branchLength, branchId), constants{constants} {

    CaResSpines.resize(branchLength / gap, CaResSynapseSpine());
    for (CaResSynapseSpine &spine : CaResSpines) {
        spine.resourcesAvailable = constants.initialResources;
        // spine.calciumFree        = constants.initialCalcium;
        spine.weight = constants.initialWeight;
        spine.PreDiffusion();
    }
    if (branchSlots < 1) {
        throw "A branch with one synapse slot cannot do diffusion";
    }
    // waitingMatrix.resize(prespikeDelaySteps, std::vector<double>(branchLength / gap, 0));
}

void CaDiffusionBranch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
    return;
}

// void CaDiffusionBranch::PreSpikeCalciumInflux(TStepInt &&step) {
//     for (size_t index : std::ranges::views::iota(0u, CaResSpines.size())) {
//         // For now, we consider that unconnected segments of the branch also have calcium influx
//         //  if (!CaResSynapseSpine.at(index).connected){
//         //    continue;
//         //  }
//         CaResSpines.at(index).calciumFree += waitingMatrix.at(step).at(index);
//     }
//     // for (size_t index : std::ranges::views::iota(0u, prespikeWaitingMatrix.size()-1)){
//     //   prespikeWaitingMatrix.at(index).swap(prespikeWaitingMatrix.at(index+1));
//     // }
//     std::fill(PAR_UNSEQ, waitingMatrix.at(step).begin(), waitingMatrix.at(step).end(), constants.calciumInfluxBasal);
// }

// void CaDiffusionBranch::Advect(TStepInt step) {
void CaDiffusionBranch::Advect() {
    // A possible ugly optimization is to do the first two iteration containers (if connected), the diffusion of the first container, then iterate from index 2 onwards and do the diffusion of
    // container index-1. When the end is reached, you do the diffusion of last container.
    double inactiveCalmodulin{}, kDot{}, nDot{}, kPDot{}, wDot{}, camDot{}, ngDot{}; // inactivePhosphatases{},inactiveKinases{}, unboundNeurogranin{}; is only used once, so you can
                                                                                     // calculate it in place
    const Constants constants{constants};                                                                                 
    size_t lastIndex{CaResSpines.size() - 1};
    // PreSpikeCalciumInflux(std::move(step));
    // I apologize in advance, as this function requires being run with spatial and temporal locality for max speed.
    // Spatial locality is achieved in the spine vector (not pointers). Temporal locality is achieved by making a mess
    // of a function
    // All constants.reactions happen in a single loop because diffusion is separate
    for (CaResSynapseSpine &spine : CaResSpines) { // If diffusion does not have to be inside here, we iterate the container
        // 1st, calcium input (other function)
        if(!spine.connected){
            continue;
        }
        //Here we influx calcium with the nonlinear traces
        spine.preTransientIncrease -=spine.preTransientIncrease*constants.preCalciumDecayRate; //This is B in Graupner and Brunel appendix
        spine.preTransient+=constants.preCalciumFluxFactor*(spine.preTransientIncrease-spine.preTransient*constants.preCalciumRiseRate);//This is A in appendix
        spine.postTransientIncrease -=spine.postTransientIncrease*constants.postCalciumDecayRate; //This is F in Graupner and Brunel appendix
        spine.postTransient+=constants.postCalciumFluxFactor*(spine.postTransientIncrease-spine.postTransient*constants.postCalciumRiseRate);//This is E in appendix
        //Here we add the traces plus the basal flux to the calcium
        spine.calciumFree+=spine.preTransient+spine.postTransient;
        //Until here, now first reactions
        inactiveCalmodulin = constants.calmodulinTotal - spine.calmodulinActive - spine.calmodulinNeurogranin;
        // 2nd, neurogranin constants.reaction (might swap for third)
        ngDot = constants.reaction1Ctt * (constants.neurograninTotal - spine.calmodulinNeurogranin) * (inactiveCalmodulin)-constants.reaction2Ctt * spine.calmodulinNeurogranin;
        spine.calmodulinNeurogranin += ngDot;
        inactiveCalmodulin -= ngDot;
        // 3rd Ca-CaM constants.reaction (might swap again) with NG consumption
        camDot = constants.reaction3Ctt * std::pow(spine.calciumFree, 2) * inactiveCalmodulin - constants.reaction4Ctt * spine.calmodulinActive;
        spine.calmodulinActive += camDot;
        spine.calciumFree -= 2 * camDot;
        // ORDERING
        //  4th calcineurin binding
        nDot = constants.reaction5Ctt * (constants.calcineurinTotal - spine.calcineurinActive) * spine.calmodulinActive - constants.reaction6Ctt * spine.calcineurinActive;
        spine.calcineurinActive += nDot;
        // 5th kinase autophosphorylation
        kPDot = constants.reaction7Ctt * spine.kinasesCaM * (spine.kinasesCaM + spine.kinasesPhospho) - constants.reaction8Ctt * spine.calcineurinActive * spine.kinasesPhospho;
        spine.kinasesPhospho += kPDot;
        spine.kinasesCaM -= kPDot;
        // 6th kinase activation via CaM -kP
        kDot = constants.reaction9Ctt * (constants.kinasesTotal - spine.kinasesCaM - spine.kinasesPhospho) - constants.reaction10Ctt * spine.kinasesCaM;
        spine.kinasesCaM += kDot;
        // ORDERING
        //  7th active CaM consumption by Ndot and Kdot (not kPdot)
        spine.calmodulinActive -= nDot + kDot;
        // 8th Change in synapse spine size/weight
        wDot = constants.reaction11Ctt * spine.resourcesAvailable * (spine.kinasesCaM + spine.kinasesPhospho) - constants.reaction12Ctt * spine.weight * spine.calcineurinActive;
        spine.weight += wDot;
        // 9th Consumption of resources by weight change
        spine.resourcesAvailable -= wDot;
        // For the next timestep
        spine.PreDiffusion();
    }
    //  Both diffusion functions are executed simultaneously for temporal locality
    // First boundary case
    // Diffusion of calcium
    CaResSpines.at(0).calciumFree += constants.caDiffusionFct * (-CaResSpines.at(0).calciumOldStep + CaResSpines.at(1).calciumOldStep);
    CaResSpines.at(0).calciumFree += constants.calciumInfluxBasal-CaResSpines.at(0).calciumFree * constants.calciumBufferingCtt;
    // Diffusion of resources
    CaResSpines.at(0).resourcesAvailable += constants.resourceDiffusionFct * (-CaResSpines.at(0).resourcesOldStep + CaResSpines.at(1).resourcesOldStep);
    for (size_t spineIndex : std::ranges::views::iota(1u, lastIndex)) {
        // Diffusion of calcium
        CaResSpines.at(spineIndex).calciumFree +=
            constants.caDiffusionFct * (-2 * CaResSpines.at(spineIndex).calciumOldStep + CaResSpines.at(spineIndex - 1).calciumOldStep + CaResSpines.at(spineIndex + 1).calciumOldStep);
        CaResSpines.at(spineIndex).calciumFree += constants.calciumInfluxBasal-CaResSpines.at(spineIndex).calciumFree * constants.calciumBufferingCtt;
        // Diffusion of resources
        CaResSpines.at(spineIndex).resourcesAvailable +=
            constants.resourceDiffusionFct * (-2 * CaResSpines.at(spineIndex).resourcesOldStep + CaResSpines.at(spineIndex - 1).resourcesOldStep + CaResSpines.at(spineIndex + 1).resourcesOldStep);
    }
    // Last boundary case
    //  Diffusion of calcium
    CaResSpines.at(lastIndex).calciumFree += constants.caDiffusionFct * (-CaResSpines.at(lastIndex).calciumOldStep + CaResSpines.at(lastIndex - 1).calciumOldStep);
    CaResSpines.at(lastIndex).calciumFree += constants.calciumInfluxBasal-CaResSpines.at(lastIndex).calciumFree * constants.calciumBufferingCtt;
    // Diffusion of resources
    CaResSpines.at(lastIndex).resourcesAvailable += constants.resourceDiffusionFct * (-CaResSpines.at(lastIndex).resourcesOldStep + CaResSpines.at(lastIndex - 1).resourcesOldStep);
}
void CaDiffusionBranch::PostSpikeCalciumFlux() {
    double nonlinearFactor{constants.nonlinearFactorNMDA};
    for (CaResSynapseSpine &spine : CaResSpines) {
        spine.postTransientIncrease+=(1.+nonlinearFactor*spine.preTransient);
    }
}
// void CaDiffusionBranch::Advect(TStepInt step) {
//     double inactiveCalmodulin{}, kDot{}, nDot{}, kPDot{}, wDot{}, camDot{},
//         ngDot{}; // inactivePhosphatases{},inactiveKinases{}, unboundNeurogranin{}; is only used once, so you can
//                  // calculate it in place
//     size_t lastIndex{CaResSpines.size()-1};
//     PreSpikeCalciumInflux(std::move(step));
//     // I apologize in advance, as this function requires being run with spatial and temporal locality for max speed.
//     // Spatial locality is achieved in the spine vector (not pointers). Temporal locality is achieved by making a mess
//     // of a function
//     //All constants.reactions happen in a single loop because diffusion is separate
//     for (size_t spineIndex : std::ranges::views::iota(0u, CaResSpines.size())) { //If diffusion does not have to be inside here, we iterate the container
//         if(!caResSpines.at(spineIndex).connected){
//         continue;
// }
// 1st, calcium input (other function)
//         inactiveCalmodulin=calmodulinTotal-CaResSpines.at(spineIndex).calmodulinActive-CaResSpines.at(spineIndex).calmodulinNeurogranin;
//         // 2nd, neurogranin constants.reaction (might swap for third)
//         ngDot = constants.reaction7Ctt*(neurograninTotal-CaResSpines.at(spineIndex).calmodulinNeurogranin)*(inactiveCalmodulin)-reaction8Ctt*CaResSpines.at(spineIndex).calmodulinNeurogranin;
//         CaResSpines.at(spineIndex).calmodulinNeurogranin+=ngDot;
//         inactiveCalmodulin-=ngDot;
//         // 3rd Ca-CaM constants.reaction (might swap again) with NG consumption
//         camDot = constants.reaction9Ctt*std::pow(CaResSpines.at(spineIndex).calciumFree,2)*inactiveCalmodulin - constants.reaction10Ctt*CaResSpines.at(spineIndex).calmodulinActive;
//         CaResSpines.at(spineIndex).calmodulinActive+=camDot;
//         CaResSpines.at(spineIndex).calciumFree-=2*camDot;
//         //ORDERING
//         // 4th calcineurin binding
//         nDot = constants.reaction3Ctt*(calcineurinTotal-CaResSpines.at(spineIndex).calcineurinActive)*CaResSpines.at(spineIndex).calmodulinActive -
//         constants.reaction4Ctt*CaResSpines.at(spineIndex).calcineurinActive; CaResSpines.at(spineIndex).calcineurinActive+=nDot;
//         // 5th kinase autophosphorylation
//         kPDot= constants.reaction11Ctt + constants.reaction12Ctt;
//         CaResSpines.at(spineIndex).kinasesPhospho+=kPDot;
//         CaResSpines.at(spineIndex).kinasesCaM-=kPDot;
//         // 6th kinase activation via CaM -kP
//         kDot= constants.reaction1Ctt + constants.reaction2Ctt;
//         CaResSpines.at(spineIndex).kinasesCaM+=kDot;
//         //ORDERING
//         // 7th active CaM consumption by Ndot and Kdot (not kPdot)
//         CaResSpines.at(spineIndex).calmodulinActive-=nDot+kDot;
//         // 8th Change in synapse spine size/weight
//         wDot=reaction5Ctt - constants.reaction6Ctt;
//         CaResSpines.at(spineIndex).weight+=wDot;
//         // 9th Consumption of resources by weight change
//         CaResSpines.at(spineIndex).resourcesAvailable-=wDot;
//         //For the next timestep
//         CaResSpines.at(spineIndex).PreDiffusion();
//     }
//     //  Both diffusion functions are executed simultaneously for temporal locality
//     //First boundary case
//           // Diffusion of calcium
//     CaResSpines.at(0).calciumFree+=caDiffusionFct*(-CaResSpines.at(0).calciumOldStep+CaResSpines.at(1).calciumOldStep);
//     CaResSpines.at(0).calciumFree-=CaResSpines.at(0).calciumFree*calciumBufferingCtt;
//           // Diffusion of resources
//     CaResSpines.at(0).resourcesAvailable+=resourceDiffusionFct*(-CaResSpines.at(0).resourcesOldStep+CaResSpines.at(1).resourcesOldStep);
//     for (size_t spineIndex : std::ranges::views::iota(1u, lastIndex)) {
//       // Diffusion of calcium
//         CaResSpines.at(spineIndex).calciumFree+=caDiffusionFct*(-2*CaResSpines.at(spineIndex).calciumOldStep+CaResSpines.at(spineIndex-1).calciumOldStep+CaResSpines.at(spineIndex+1).calciumOldStep);
//         CaResSpines.at(spineIndex).calciumFree-=CaResSpines.at(spineIndex).calciumFree*calciumBufferingCtt;
//               // Diffusion of resources
//         CaResSpines.at(spineIndex).resourcesAvailable+=resourceDiffusionFct*(-2*CaResSpines.at(spineIndex).resourcesOldStep+CaResSpines.at(spineIndex-1).resourcesOldStep+CaResSpines.at(spineIndex+1).resourcesOldStep);
//     }
//     //Last boundary case
//               // Diffusion of calcium
//     CaResSpines.at(lastIndex).calciumFree+=caDiffusionFct*(-CaResSpines.at(lastIndex).calciumOldStep+CaResSpines.at(lastIndex-1).calciumOldStep);
//     CaResSpines.at(lastIndex).calciumFree-=CaResSpines.at(lastIndex).calciumFree*calciumBufferingCtt;
//                   // Diffusion of resources
//     CaResSpines.at(lastIndex).resourcesAvailable+=resourceDiffusionFct*(-CaResSpines.at(lastIndex).resourcesOldStep+CaResSpines.at(lastIndex-1).resourcesOldStep);
// }
