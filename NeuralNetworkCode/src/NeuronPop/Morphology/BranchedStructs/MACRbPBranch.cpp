#include "./MACRbpBranch.hpp"
#include "MACRbPBranch.hpp"

MACRbPBranch::MACRbPBranch(std::vector<int> anteriorBranches, double gap, double branchLength, int branchId, Constants constants,
                           MACRbPSynapseSpine spine)
    : Branch(anteriorBranches, gap, branchLength, branchId), constants{constants} {
  MACRbPspines.resize(branchLength / gap, spine);
  // for (MACRbPSynapseSpine &spine : CaResSpines) {
  //   spine.resourcesAvailable = constants.initialResources;
  //   spine.weight             = constants.initialWeight;
  // }
  // waitingMatrix.resize(prespikeDelaySteps, std::vector<double>(branchLength / gap, constants.calciumInfluxBasal));
}

MACRbPBranch::MACRbPBranch(double gap, double branchLength, int branchId, Constants constants, MACRbPSynapseSpine spine)
    : Branch(gap, branchLength, branchId), constants{constants} {

  MACRbPspines.resize(branchLength / gap, spine);
  for (MACRbPSynapseSpine &spine : MACRbPspines) {
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
  // If it is ever necessary to release the weight of unconnected synapses as free resources, uncomment the next lines of code.
  // PROBLEM WITH THIS: Spines will not be in steady state anymore, they will start over-depressed
  // double extraResources{};
  // for (const MACRbPSynapseSpine &spine : MACRbPspines) {
  //   if (!spine.connected) {
  //     extraResources += spine.weight / constants.resourceConversionFct;
  //   }
  // }
  // double extra_alloc_resources{extraResources / static_cast<double>(MACRbPspines.size())};
  // // for (MACRbPSynapseSpine &spine : MACRbPspines) {
  // //   spine.resourcesAvailable += extra_alloc_resources;
  // // }
  // std::transform(PAR_UNSEQ, MACRbPspines.begin(), MACRbPspines.end(), MACRbPspines.begin(),
  //                [extra_alloc_resources](MACRbPSynapseSpine &spine) { spine.resourcesAvailable += extra_alloc_resources; });
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
  std::for_each(PAR_UNSEQ, MACRbPspines.begin(), MACRbPspines.end(), [ctt](MACRbPSynapseSpine &spine) {
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
    }
    spine.PreDiffusion();
#ifndef NDEBUG
    spine.CheckNegativeValues(ctt);
#endif
  });
  // DIFFUSION STARTS HERE
  MACRbPspines.at(0).calciumFree += ctt.caDiffusionFct * (-MACRbPspines.at(0).calciumOldStep + MACRbPspines.at(1).calciumOldStep);
  MACRbPspines.at(0).calciumFree += ctt.calciumInfluxBasal - MACRbPspines.at(0).calciumFree * ctt.calciumExtrusionCtt;
  // Diffusion of resources
  MACRbPspines.at(0).resourcesAvailable += ctt.resourceDiffusionFct * (-MACRbPspines.at(0).resourcesOldStep + MACRbPspines.at(1).resourcesOldStep);

  size_t lastIndex{MACRbPspines.size() - 1};
  auto   range = std::ranges::common_view(std::views::iota(1ull) | std::views::take(lastIndex));
  std::for_each(PAR_UNSEQ, range.begin(), range.end(), [ctt, this](size_t spineIndex) {
    // Diffusion of calcium
    MACRbPspines.at(spineIndex).calciumFree +=
        ctt.caDiffusionFct * (-2 * MACRbPspines.at(spineIndex).calciumOldStep + MACRbPspines.at(spineIndex - 1).calciumOldStep +
                              MACRbPspines.at(spineIndex + 1).calciumOldStep);
    MACRbPspines.at(spineIndex).calciumFree += ctt.calciumInfluxBasal - MACRbPspines.at(spineIndex).calciumFree * ctt.calciumExtrusionCtt;
    // Diffusion of resources
    MACRbPspines.at(spineIndex).resourcesAvailable +=
        ctt.resourceDiffusionFct * (-2 * MACRbPspines.at(spineIndex).resourcesOldStep + MACRbPspines.at(spineIndex - 1).resourcesOldStep +
                                    MACRbPspines.at(spineIndex + 1).resourcesOldStep);
  });
  // Last boundary case
  //  Diffusion of calcium
  MACRbPspines.at(lastIndex).calciumFree +=
      ctt.caDiffusionFct * (-MACRbPspines.at(lastIndex).calciumOldStep + MACRbPspines.at(lastIndex - 1).calciumOldStep);
  MACRbPspines.at(lastIndex).calciumFree += ctt.calciumInfluxBasal - MACRbPspines.at(lastIndex).calciumFree * ctt.calciumExtrusionCtt;
  // Diffusion of resources
  MACRbPspines.at(lastIndex).resourcesAvailable +=
      ctt.resourceDiffusionFct * (-MACRbPspines.at(lastIndex).resourcesOldStep + MACRbPspines.at(lastIndex - 1).resourcesOldStep);
}

void MACRbPBranch::AdvectUnrolled() {
  const Constants ctt{this->constants};
  // Calc+resources STARTS HERE
  MACRbPspines.at(0).PreDiffusion(), MACRbPspines.at(1).PreDiffusion();
  MACRbPspines.at(0).calciumFree += ctt.caDiffusionFct * (-MACRbPspines.at(0).calciumOldStep + MACRbPspines.at(1).calciumOldStep);
  MACRbPspines.at(0).calciumFree += ctt.calciumInfluxBasal - MACRbPspines.at(0).calciumFree * ctt.calciumExtrusionCtt;
  // Diffusion of resources
  MACRbPspines.at(0).resourcesAvailable += ctt.resourceDiffusionFct * (-MACRbPspines.at(0).resourcesOldStep + MACRbPspines.at(1).resourcesOldStep);
  // Reactions
  if (MACRbPspines.at(0).connected) {
    // Here we influx calcium with the nonlinear traces
    MACRbPspines.at(0).preTransientIncrease -=
        MACRbPspines.at(0).preTransientIncrease * ctt.preCalciumDecayRate; // This is B in Graupner and Brunel appendix
    MACRbPspines.at(0).preTransient += ctt.preCalciumFluxFactor * (MACRbPspines.at(0).preTransientIncrease -
                                                                   MACRbPspines.at(0).preTransient * ctt.preCalciumRiseRate); // This is A in appendix

    MACRbPspines.at(0).postTransientIncrease -=
        MACRbPspines.at(0).postTransientIncrease * ctt.postCalciumDecayRate; // This is F in Graupner and Brunel appendix
    MACRbPspines.at(0).postTransient +=
        ctt.postCalciumFluxFactor *
        (MACRbPspines.at(0).postTransientIncrease - MACRbPspines.at(0).postTransient * ctt.postCalciumRiseRate); // This is E in appendix
    // Here we add the traces plus the basal flux to the calcium
    MACRbPspines.at(0).calciumFree += MACRbPspines.at(0).preTransient + MACRbPspines.at(0).postTransient;
    // Until here, now first reactions

    double freeNg{ctt.neurograninTotal - MACRbPspines.at(0).calmodulinNeurogranin}, epsilon{};
    for (int i = 0; i < ctt.newtonIterations; i++) {
      epsilon = (ctt.reaction1234Ctt * freeNg * MACRbPspines.at(0).calmodulinActive -
                 MACRbPspines.at(0).calmodulinNeurogranin * std::pow(MACRbPspines.at(0).calciumFree, 2)) /
                (ctt.reaction1234Ctt * (freeNg + MACRbPspines.at(0).calmodulinActive) +
                 4 * MACRbPspines.at(0).calmodulinNeurogranin * MACRbPspines.at(0).calciumFree + std::pow(MACRbPspines.at(0).calciumFree, 2));
      freeNg -= epsilon;
      MACRbPspines.at(0).calmodulinActive -= epsilon;
      MACRbPspines.at(0).calmodulinNeurogranin += epsilon;
      MACRbPspines.at(0).calciumFree += 2 * epsilon;
    }

    double kDot{}, nDot{}, kPDot{}, wDot{};
    nDot = ctt.reaction5Ctt * (ctt.calcineurinTotal - MACRbPspines.at(0).calcineurinActive) * MACRbPspines.at(0).calmodulinActive -
           ctt.reaction6Ctt * MACRbPspines.at(0).calcineurinActive;
    MACRbPspines.at(0).calcineurinActive += nDot;
    MACRbPspines.at(0).calmodulinActive -= nDot;
    // 5th kinase autophosphorylation
    kPDot = ctt.reaction7Ctt * MACRbPspines.at(0).kinasesCaM * (MACRbPspines.at(0).kinasesCaM + MACRbPspines.at(0).kinasesPhospho) -
            ctt.reaction8Ctt * MACRbPspines.at(0).calcineurinActive * MACRbPspines.at(0).kinasesPhospho;
    MACRbPspines.at(0).kinasesPhospho += kPDot;
    MACRbPspines.at(0).kinasesCaM -= kPDot;
    // 6th kinase activation via CaM -kP
    kDot = ctt.reaction9Ctt * (ctt.kinasesTotal - MACRbPspines.at(0).kinasesCaM - MACRbPspines.at(0).kinasesPhospho) *
               MACRbPspines.at(0).calmodulinActive -
           ctt.reaction10Ctt * MACRbPspines.at(0).kinasesCaM;
    MACRbPspines.at(0).kinasesCaM += kDot;
    MACRbPspines.at(0).calmodulinActive -= kDot;
    // ORDERING
    //  7th active CaM consumption by Ndot and Kdot (not kPdot)
    // spine.calmodulinActive -= nDot + kDot;//We are not doing this because this could mean negative active calmodulins
    // 8th Change in synapse spine size/weight
    wDot = ctt.reaction11Ctt * MACRbPspines.at(0).resourcesAvailable * (MACRbPspines.at(0).kinasesCaM + MACRbPspines.at(0).kinasesPhospho) -
           ctt.reaction12Ctt * MACRbPspines.at(0).weight * MACRbPspines.at(0).calcineurinActive;
    MACRbPspines.at(0).weight += wDot;
    // 9th Consumption of resources by weight change
    MACRbPspines.at(0).resourcesAvailable -=
        (wDot) / ctt.resourceConversionFct; // This should be the case for converting the dmV/spike to concentration of resources
    // For the next timestep
  }
#ifndef NDEBUG
  MACRbPspines.at(0).CheckNegativeValues(ctt);
#endif
  size_t lastIndex{MACRbPspines.size() - 1};
  auto   range = std::ranges::common_view(std::views::iota(1ull) | std::views::take(lastIndex)); // It takes up to that element
  std::for_each(PAR, range.begin(), range.end(), [this, ctt](size_t spineIndex) {
    MACRbPspines.at(spineIndex + 1).PreDiffusion();
    // Diffusion of calcium
    MACRbPspines.at(spineIndex).calciumFree +=
        ctt.caDiffusionFct * (-2 * MACRbPspines.at(spineIndex).calciumOldStep + MACRbPspines.at(spineIndex - 1).calciumOldStep +
                              MACRbPspines.at(spineIndex + 1).calciumOldStep);
    MACRbPspines.at(spineIndex).calciumFree += ctt.calciumInfluxBasal - MACRbPspines.at(spineIndex).calciumFree * ctt.calciumExtrusionCtt;
    // Diffusion of resources
    MACRbPspines.at(spineIndex).resourcesAvailable +=
        ctt.resourceDiffusionFct * (-2 * MACRbPspines.at(spineIndex).resourcesOldStep + MACRbPspines.at(spineIndex - 1).resourcesOldStep +
                                    MACRbPspines.at(spineIndex + 1).resourcesOldStep);
    if (MACRbPspines.at(spineIndex).connected) {
      // Here we influx calcium with the nonlinear traces
      MACRbPspines.at(spineIndex).preTransientIncrease -=
          MACRbPspines.at(spineIndex).preTransientIncrease * ctt.preCalciumDecayRate; // This is B in Graupner and Brunel appendix
      MACRbPspines.at(spineIndex).preTransient +=
          ctt.preCalciumFluxFactor * (MACRbPspines.at(spineIndex).preTransientIncrease -
                                      MACRbPspines.at(spineIndex).preTransient * ctt.preCalciumRiseRate); // This is A in appendix

      MACRbPspines.at(spineIndex).postTransientIncrease -=
          MACRbPspines.at(spineIndex).postTransientIncrease * ctt.postCalciumDecayRate; // This is F in Graupner and Brunel appendix
      MACRbPspines.at(spineIndex).postTransient +=
          ctt.postCalciumFluxFactor * (MACRbPspines.at(spineIndex).postTransientIncrease -
                                       MACRbPspines.at(spineIndex).postTransient * ctt.postCalciumRiseRate); // This is E in appendix
      // Here we add the traces plus the basal flux to the calcium
      MACRbPspines.at(spineIndex).calciumFree += MACRbPspines.at(spineIndex).preTransient + MACRbPspines.at(spineIndex).postTransient;
      // Until here, now first reactions

      double freeNg{ctt.neurograninTotal - MACRbPspines.at(spineIndex).calmodulinNeurogranin}, epsilon{};
      for (int i = 0; i < ctt.newtonIterations; i++) {
        epsilon = (ctt.reaction1234Ctt * freeNg * MACRbPspines.at(spineIndex).calmodulinActive -
                   MACRbPspines.at(spineIndex).calmodulinNeurogranin * std::pow(MACRbPspines.at(spineIndex).calciumFree, 2)) /
                  (ctt.reaction1234Ctt * (freeNg + MACRbPspines.at(spineIndex).calmodulinActive) +
                   4 * MACRbPspines.at(spineIndex).calmodulinNeurogranin * MACRbPspines.at(spineIndex).calciumFree +
                   std::pow(MACRbPspines.at(spineIndex).calciumFree, 2));
        freeNg -= epsilon;
        MACRbPspines.at(spineIndex).calmodulinActive -= epsilon;
        MACRbPspines.at(spineIndex).calmodulinNeurogranin += epsilon;
        MACRbPspines.at(spineIndex).calciumFree += 2 * epsilon;
      }

      double kDot{}, nDot{}, kPDot{}, wDot{};
      nDot =
          ctt.reaction5Ctt * (ctt.calcineurinTotal - MACRbPspines.at(spineIndex).calcineurinActive) * MACRbPspines.at(spineIndex).calmodulinActive -
          ctt.reaction6Ctt * MACRbPspines.at(spineIndex).calcineurinActive;
      MACRbPspines.at(spineIndex).calcineurinActive += nDot;
      MACRbPspines.at(spineIndex).calmodulinActive -= nDot;
      // 5th kinase autophosphorylation
      kPDot = ctt.reaction7Ctt * MACRbPspines.at(spineIndex).kinasesCaM *
                  (MACRbPspines.at(spineIndex).kinasesCaM + MACRbPspines.at(spineIndex).kinasesPhospho) -
              ctt.reaction8Ctt * MACRbPspines.at(spineIndex).calcineurinActive * MACRbPspines.at(spineIndex).kinasesPhospho;
      MACRbPspines.at(spineIndex).kinasesPhospho += kPDot;
      MACRbPspines.at(spineIndex).kinasesCaM -= kPDot;
      // 6th kinase activation via CaM -kP
      kDot = ctt.reaction9Ctt * (ctt.kinasesTotal - MACRbPspines.at(spineIndex).kinasesCaM - MACRbPspines.at(spineIndex).kinasesPhospho) *
                 MACRbPspines.at(spineIndex).calmodulinActive -
             ctt.reaction10Ctt * MACRbPspines.at(spineIndex).kinasesCaM;
      MACRbPspines.at(spineIndex).kinasesCaM += kDot;
      MACRbPspines.at(spineIndex).calmodulinActive -= kDot;
      // ORDERING
      //  7th active CaM consumption by Ndot and Kdot (not kPdot)
      // spine.calmodulinActive -= nDot + kDot;//We are not doing this because this could mean negative active calmodulins
      // 8th Change in synapse spine size/weight
      wDot = ctt.reaction11Ctt * MACRbPspines.at(spineIndex).resourcesAvailable *
                 (MACRbPspines.at(spineIndex).kinasesCaM + MACRbPspines.at(spineIndex).kinasesPhospho) -
             ctt.reaction12Ctt * MACRbPspines.at(spineIndex).weight * MACRbPspines.at(spineIndex).calcineurinActive;
      MACRbPspines.at(spineIndex).weight += wDot;
      // 9th Consumption of resources by weight change
      MACRbPspines.at(spineIndex).resourcesAvailable -=
          (wDot) / ctt.resourceConversionFct; // This should be the case for converting the dmV/spike to concentration of resources
      // For the next timestep
    }
#ifndef NDEBUG
    MACRbPspines.at(spineIndex).CheckNegativeValues(ctt);
#endif
  });
  // Last boundary case
  //  Diffusion of calcium
  MACRbPspines.at(lastIndex).calciumFree +=
      ctt.caDiffusionFct * (-MACRbPspines.at(lastIndex).calciumOldStep + MACRbPspines.at(lastIndex - 1).calciumOldStep);
  MACRbPspines.at(lastIndex).calciumFree += ctt.calciumInfluxBasal - MACRbPspines.at(lastIndex).calciumFree * ctt.calciumExtrusionCtt;
  // Diffusion of resources
  MACRbPspines.at(lastIndex).resourcesAvailable +=
      ctt.resourceDiffusionFct * (-MACRbPspines.at(lastIndex).resourcesOldStep + MACRbPspines.at(lastIndex - 1).resourcesOldStep);
  // Reactions
  if (MACRbPspines.at(lastIndex).connected) {
    // Here we influx calcium with the nonlinear traces
    MACRbPspines.at(lastIndex).preTransientIncrease -=
        MACRbPspines.at(lastIndex).preTransientIncrease * ctt.preCalciumDecayRate; // This is B in Graupner and Brunel appendix
    MACRbPspines.at(lastIndex).preTransient +=
        ctt.preCalciumFluxFactor *
        (MACRbPspines.at(lastIndex).preTransientIncrease - MACRbPspines.at(lastIndex).preTransient * ctt.preCalciumRiseRate); // This is A in appendix

    MACRbPspines.at(lastIndex).postTransientIncrease -=
        MACRbPspines.at(lastIndex).postTransientIncrease * ctt.postCalciumDecayRate; // This is F in Graupner and Brunel appendix
    MACRbPspines.at(lastIndex).postTransient +=
        ctt.postCalciumFluxFactor * (MACRbPspines.at(lastIndex).postTransientIncrease -
                                     MACRbPspines.at(lastIndex).postTransient * ctt.postCalciumRiseRate); // This is E in appendix
    // Here we add the traces plus the basal flux to the calcium
    MACRbPspines.at(lastIndex).calciumFree += MACRbPspines.at(lastIndex).preTransient + MACRbPspines.at(lastIndex).postTransient;
    // Until here, now first reactions

    double freeNg{ctt.neurograninTotal - MACRbPspines.at(lastIndex).calmodulinNeurogranin}, epsilon{};
    for (int i = 0; i < ctt.newtonIterations; i++) {
      epsilon = (ctt.reaction1234Ctt * freeNg * MACRbPspines.at(lastIndex).calmodulinActive -
                 MACRbPspines.at(lastIndex).calmodulinNeurogranin * std::pow(MACRbPspines.at(lastIndex).calciumFree, 2)) /
                (ctt.reaction1234Ctt * (freeNg + MACRbPspines.at(lastIndex).calmodulinActive) +
                 4 * MACRbPspines.at(lastIndex).calmodulinNeurogranin * MACRbPspines.at(lastIndex).calciumFree +
                 std::pow(MACRbPspines.at(lastIndex).calciumFree, 2));
      freeNg -= epsilon;
      MACRbPspines.at(lastIndex).calmodulinActive -= epsilon;
      MACRbPspines.at(lastIndex).calmodulinNeurogranin += epsilon;
      MACRbPspines.at(lastIndex).calciumFree += 2 * epsilon;
    }

    double kDot{}, nDot{}, kPDot{}, wDot{};
    nDot = ctt.reaction5Ctt * (ctt.calcineurinTotal - MACRbPspines.at(lastIndex).calcineurinActive) * MACRbPspines.at(lastIndex).calmodulinActive -
           ctt.reaction6Ctt * MACRbPspines.at(lastIndex).calcineurinActive;
    MACRbPspines.at(lastIndex).calcineurinActive += nDot;
    MACRbPspines.at(lastIndex).calmodulinActive -= nDot;
    // 5th kinase autophosphorylation
    kPDot = ctt.reaction7Ctt * MACRbPspines.at(lastIndex).kinasesCaM *
                (MACRbPspines.at(lastIndex).kinasesCaM + MACRbPspines.at(lastIndex).kinasesPhospho) -
            ctt.reaction8Ctt * MACRbPspines.at(lastIndex).calcineurinActive * MACRbPspines.at(lastIndex).kinasesPhospho;
    MACRbPspines.at(lastIndex).kinasesPhospho += kPDot;
    MACRbPspines.at(lastIndex).kinasesCaM -= kPDot;
    // 6th kinase activation via CaM -kP
    kDot = ctt.reaction9Ctt * (ctt.kinasesTotal - MACRbPspines.at(lastIndex).kinasesCaM - MACRbPspines.at(lastIndex).kinasesPhospho) *
               MACRbPspines.at(lastIndex).calmodulinActive -
           ctt.reaction10Ctt * MACRbPspines.at(lastIndex).kinasesCaM;
    MACRbPspines.at(lastIndex).kinasesCaM += kDot;
    MACRbPspines.at(lastIndex).calmodulinActive -= kDot;
    // ORDERING
    //  7th active CaM consumption by Ndot and Kdot (not kPdot)
    // spine.calmodulinActive -= nDot + kDot;//We are not doing this because this could mean negative active calmodulins
    // 8th Change in synapse spine size/weight
    wDot = ctt.reaction11Ctt * MACRbPspines.at(lastIndex).resourcesAvailable *
               (MACRbPspines.at(lastIndex).kinasesCaM + MACRbPspines.at(lastIndex).kinasesPhospho) -
           ctt.reaction12Ctt * MACRbPspines.at(lastIndex).weight * MACRbPspines.at(lastIndex).calcineurinActive;
    MACRbPspines.at(lastIndex).weight += wDot;
    // 9th Consumption of resources by weight change
    MACRbPspines.at(lastIndex).resourcesAvailable -=
        (wDot) / ctt.resourceConversionFct; // This should be the case for converting the dmV/spike to concentration of resources
    // For the next timestep
  }
#ifndef NDEBUG
  MACRbPspines.at(lastIndex).CheckNegativeValues(ctt);
#endif
}

double MACRbPBranch::GetTotalWeight() const {
  return std::reduce(MACRbPspines.begin(), MACRbPspines.end(), 0.0, [](double accumulator, const MACRbPSynapseSpine &spine) {
    if (spine.connected)
      return accumulator + spine.weight;
    else
      return accumulator;
  });
}
void MACRbPBranch::PostSpikeCalciumFlux() {
  const double nonlinearFactor{constants.nonlinearFactorNMDA};
  for (MACRbPSynapseSpine &spine : MACRbPspines) {
    spine.postTransientIncrease += (1. + nonlinearFactor * spine.preTransient);
  }
}