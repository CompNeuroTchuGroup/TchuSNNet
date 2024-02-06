#include "./MACRbpBranch.hpp"

MACRbPBranch::MACRbPBranch(std::vector<int> anteriorBranches, double gap, double branchLength, int branchId, Constants constants)
    : Branch(anteriorBranches, gap, branchLength, branchId), constants{constants} {
  CaResSpines.resize(branchLength / gap, MACRbPSynapseSpine());
  for (MACRbPSynapseSpine &spine : CaResSpines) {
    spine.resourcesAvailable = constants.initialResources;
    spine.weight             = constants.initialWeight;
  }
  // waitingMatrix.resize(prespikeDelaySteps, std::vector<double>(branchLength / gap, constants.calciumInfluxBasal));
}

MACRbPBranch::MACRbPBranch(double gap, double branchLength, int branchId, Constants constants)
    : Branch(gap, branchLength, branchId), constants{constants} {

  CaResSpines.resize(branchLength / gap, MACRbPSynapseSpine());
  for (MACRbPSynapseSpine &spine : CaResSpines) {
    spine.resourcesAvailable = constants.initialResources;
    spine.weight             = constants.initialWeight;
    spine.PreDiffusion();
  }
  if (branchSlots < 1) {
    throw "A branch with one synapse slot cannot do diffusion";
  }
  // waitingMatrix.resize(prespikeDelaySteps, std::vector<double>(branchLength / gap, 0));
}

void MACRbPBranch::PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) {
  return;
}

void MACRbPBranch::Advect() {
  // A possible ugly optimization is to do the first two iteration containers (if connected), the diffusion of the first container, then iterate
  // from index 2 onwards and do the diffusion of container index-1. When the end is reached, you do the diffusion of last container.
  const Constants ctt{this->constants}; //
  // PreSpikeCalciumInflux(std::move(step));
  // I apologize in advance, as this function requires being run with spatial and temporal locality for max speed.
  // Spatial locality is achieved in the spine vector (not pointers). Temporal locality is achieved by making a mess
  // of a function
  // All constants.reactions happen in a single loop because diffusion is separate
  std::for_each(PAR, CaResSpines.begin(), CaResSpines.end(), [ctt](MACRbPSynapseSpine &spine) {
    // const Constants ctt{ctt};
    if (spine.connected) {
      // Here we influx calcium with the nonlinear traces
      spine.preTransientIncrease -= spine.preTransientIncrease * ctt.preCalciumDecayRate; // This is B in Graupner and Brunel appendix
      spine.preTransient +=
          ctt.preCalciumFluxFactor * (spine.preTransientIncrease - spine.preTransient * ctt.preCalciumRiseRate); // This is A in appendix

      spine.postTransientIncrease -= spine.postTransientIncrease * ctt.postCalciumDecayRate; // This is F in Graupner and Brunel appendix
      spine.postTransient +=
          ctt.postCalciumFluxFactor * (spine.postTransientIncrease - spine.postTransient * ctt.postCalciumRiseRate); // This is E in appendix
      // Here we add the traces plus the basal flux to the calcium
      spine.calciumFree += spine.preTransient + spine.postTransient;
      // Until here, now first reactions

      double freeNg{ctt.neurograninTotal - spine.calmodulinNeurogranin}, epsilon{};
      for (int i = 0; i < ctt.newtonIterations; i++) {
        epsilon = (ctt.reaction1234Ctt * freeNg * spine.calmodulinActive - spine.calmodulinNeurogranin * std::pow(spine.calciumFree, 2)) /
                  (ctt.reaction1234Ctt * (freeNg + spine.calmodulinActive) + 4 * spine.calmodulinNeurogranin * spine.calciumFree +
                   std::pow(spine.calciumFree, 2));
        freeNg -= epsilon;
        spine.calmodulinActive -= epsilon;
        spine.calmodulinNeurogranin += epsilon;
        spine.calciumFree += 2 * epsilon;
      }

      // // From here
      // inactiveCalmodulin = ctt.calmodulinTotal - spine.calmodulinActive - spine.calmodulinNeurogranin;
      // //  2nd, neurogranin constants.reaction (might swap for third)
      // ngDot = ctt.reaction1Ctt * (ctt.neurograninTotal - spine.calmodulinNeurogranin) * (inactiveCalmodulin)-ctt.reaction2Ctt *
      //         spine.calmodulinNeurogranin;
      // spine.calmodulinNeurogranin += ngDot;
      // inactiveCalmodulin -= ngDot;
      // // 3rd Ca-CaM constants.reaction (might swap again) with NG consumption
      // camDot = ctt.reaction3Ctt * std::pow(spine.calciumFree, 2) * inactiveCalmodulin - ctt.reaction4Ctt * spine.calmodulinActive;
      // spine.calmodulinActive += camDot;
      // spine.calciumFree -= 2 * camDot;
      // // To here is the equilibrium
      //  ORDERING
      //   4th calcineurin binding
      double kDot{}, nDot{}, kPDot{}, wDot{};
      nDot =
          ctt.reaction5Ctt * (ctt.calcineurinTotal - spine.calcineurinActive) * spine.calmodulinActive - ctt.reaction6Ctt * spine.calcineurinActive;
      spine.calcineurinActive += nDot;
      spine.calmodulinActive -= nDot;
      // 5th kinase autophosphorylation
      kPDot = ctt.reaction7Ctt * spine.kinasesCaM * (spine.kinasesCaM + spine.kinasesPhospho) -
              ctt.reaction8Ctt * spine.calcineurinActive * spine.kinasesPhospho;
      spine.kinasesPhospho += kPDot;
      spine.kinasesCaM -= kPDot;
      // 6th kinase activation via CaM -kP
      kDot = ctt.reaction9Ctt * (ctt.kinasesTotal - spine.kinasesCaM - spine.kinasesPhospho) * spine.calmodulinActive -
             ctt.reaction10Ctt * spine.kinasesCaM;
      spine.kinasesCaM += kDot;
      spine.calmodulinActive -= kDot;
      // ORDERING
      //  7th active CaM consumption by Ndot and Kdot (not kPdot)
      // spine.calmodulinActive -= nDot + kDot;//We are not doing this because this could mean negative active calmodulins
      // 8th Change in synapse spine size/weight
      wDot = ctt.reaction11Ctt * spine.resourcesAvailable * (spine.kinasesCaM + spine.kinasesPhospho) -
             ctt.reaction12Ctt * spine.weight * spine.calcineurinActive;
      spine.weight += wDot;
      // 9th Consumption of resources by weight change
      spine.resourcesAvailable -=
          (wDot) / ctt.resourceConversionFct; // This should be the case for converting the dmV/spike to concentration of resources
      // For the next timestep
      spine.PreDiffusion();
#ifndef NDEBUG
      spine.CheckNegativeValues(ctt);
#endif
    }
  });
  // DIFFUSION STARTS HERE
  CaResSpines.at(0).calciumFree += ctt.caDiffusionFct * (-CaResSpines.at(0).calciumOldStep + CaResSpines.at(1).calciumOldStep);
  CaResSpines.at(0).calciumFree += ctt.calciumInfluxBasal - CaResSpines.at(0).calciumFree * ctt.calciumExtrusionCtt;
  // Diffusion of resources
  CaResSpines.at(0).resourcesAvailable += ctt.resourceDiffusionFct * (-CaResSpines.at(0).resourcesOldStep + CaResSpines.at(1).resourcesOldStep);

  size_t lastIndex{CaResSpines.size() - 1};
  auto   range = std::ranges::common_view(std::views::iota(1ull) | std::views::take(lastIndex - 1));
  std::for_each(PAR, range.begin(), range.end(), [this](size_t spineIndex) {
    const Constants CONST{this->constants};
    // Diffusion of calcium
    CaResSpines.at(spineIndex).calciumFree +=
        CONST.caDiffusionFct * (-2 * CaResSpines.at(spineIndex).calciumOldStep + CaResSpines.at(spineIndex - 1).calciumOldStep +
                                CaResSpines.at(spineIndex + 1).calciumOldStep);
    CaResSpines.at(spineIndex).calciumFree += CONST.calciumInfluxBasal - CaResSpines.at(spineIndex).calciumFree * CONST.calciumExtrusionCtt;
    // Diffusion of resources
    CaResSpines.at(spineIndex).resourcesAvailable +=
        CONST.resourceDiffusionFct * (-2 * CaResSpines.at(spineIndex).resourcesOldStep + CaResSpines.at(spineIndex - 1).resourcesOldStep +
                                      CaResSpines.at(spineIndex + 1).resourcesOldStep);
  });
  // Last boundary case
  //  Diffusion of calcium
  CaResSpines.at(lastIndex).calciumFree +=
      ctt.caDiffusionFct * (-CaResSpines.at(lastIndex).calciumOldStep + CaResSpines.at(lastIndex - 1).calciumOldStep);
  CaResSpines.at(lastIndex).calciumFree += ctt.calciumInfluxBasal - CaResSpines.at(lastIndex).calciumFree * ctt.calciumExtrusionCtt;
  // Diffusion of resources
  CaResSpines.at(lastIndex).resourcesAvailable +=
      ctt.resourceDiffusionFct * (-CaResSpines.at(lastIndex).resourcesOldStep + CaResSpines.at(lastIndex - 1).resourcesOldStep);
}
void MACRbPBranch::PostSpikeCalciumFlux() {
  const double nonlinearFactor{constants.nonlinearFactorNMDA};
  for (MACRbPSynapseSpine &spine : CaResSpines) {
    spine.postTransientIncrease += (1. + nonlinearFactor * spine.preTransient);
  }
}