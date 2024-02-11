#include "./MACRbPModel.hpp"
#include "MACRbPModel.hpp"

MACRbPModel::MACRbPModel(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters, MACRbPSynapseSpine &spine)
    : BranchedMorphology(infoGlobal, morphologyParameters) {
  if (spine.calmodulinActive == 0) {
    ComputeSteadyState(spine);
  }
  steady_state_spine         = spine;
  constants.newtonIterations = 2;
}

void MACRbPModel::LoadParameters(const std::vector<FileEntry> &morphologyParameters) {
  BranchedMorphology::LoadParameters(morphologyParameters);
  for (auto &[parameterName, parameterValues] : morphologyParameters) {
    if (parameterName.find("kinasesTotal") != std::string::npos) {
      this->constants.kinasesTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("phosphatasesTotal") != std::string::npos) {
      this->constants.calcineurinTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("calmodulinTotal") != std::string::npos) {
      this->constants.calmodulinTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("neurograninTotal") != std::string::npos) {
      this->constants.neurograninTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("kOne") != std::string::npos) {
      this->constants.reaction1Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kTwo") != std::string::npos) {
      this->constants.reaction2Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kThree") != std::string::npos) {
      this->constants.reaction3Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kFour") != std::string::npos) {
      this->constants.reaction4Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kFive") != std::string::npos) {
      this->constants.reaction5Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kSix") != std::string::npos) {
      this->constants.reaction6Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kSeven") != std::string::npos) {
      this->constants.reaction7Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kEight") != std::string::npos) {
      this->constants.reaction8Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kNine") != std::string::npos) {
      this->constants.reaction9Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kTen") != std::string::npos) {
      this->constants.reaction10Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kEleven") != std::string::npos) {
      this->constants.reaction11Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kTwelve") != std::string::npos) {
      this->constants.reaction12Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("resourceDiffusion") != std::string::npos) {
      this->constants.resourceDiffusionFct = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2);
    } else if (parameterName.find("calciumDiffusion") != std::string::npos) {
      this->constants.caDiffusionFct = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2);
    } else if (parameterName.find("calciumExtrusion") != std::string::npos) {
      this->constants.calciumExtrusionCtt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("initialWeight") != std::string::npos) {
      this->constants.initialWeight = std::stod(parameterValues.at(0));
    } else if (parameterName.find("initialResources") != std::string::npos) {
      this->availResourcesRatio = std::stod(parameterValues.at(0));
    } else if (parameterName.find("prespikeCalcium") != std::string::npos) {
      this->prespikeCalcium = std::stod(parameterValues.at(0));
    } else if (parameterName.find("postspikeCalcium") != std::string::npos) {
      this->postspikeCalcium = std::stod(parameterValues.at(0));
    } else if (parameterName.find("preCalciumRiseTau") != std::string::npos) {
      this->preCalciumRiseTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("preCalciumDecayTau") != std::string::npos) {
      this->preCalciumDecayTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("postCalciumRiseTau") != std::string::npos) {
      this->postCalciumRiseTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("postCalciumDecayTau") != std::string::npos) {
      this->postCalciumDecayTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("nonlinearNMDAFactor") != std::string::npos) {
      this->constants.nonlinearFactorNMDA = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("calciumBasal") != std::string::npos) {
      this->calciumBasal = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("resourceConversionFct") != std::string::npos) {
      this->constants.resourceConversionFct = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("newtonIterations") != std::string::npos) {
      this->constants.newtonIterations = (std::stoi(parameterValues.at(0)));
    }
  }
  PreCalculateParameters();
}

void MACRbPModel::PreCalculateParameters() {
  this->constants.reaction1234Ctt =
      this->constants.reaction1Ctt * this->constants.reaction4Ctt / (this->constants.reaction2Ctt * this->constants.reaction3Ctt);
  this->constants.calciumInfluxBasal = calciumBasal * constants.calciumExtrusionCtt;
  this->constants.initialResources   = this->constants.initialWeight / constants.resourceConversionFct / availResourcesRatio; // REVIEW
  this->constants.reaction11Ctt *= constants.resourceConversionFct;                                                           // REVIEW
  this->constants.preCalciumRiseRate   = 1 / preCalciumRiseTau;
  this->constants.preCalciumDecayRate  = 1 / preCalciumDecayTau;
  this->constants.postCalciumRiseRate  = 1 / postCalciumRiseTau;
  this->constants.postCalciumDecayRate = 1 / postCalciumDecayTau;

  auto amplitudeNormalization = [](double amplitude, double calcium_tau, double decay_tau) { // The taus need to be in the proper time unit
    return amplitude / (calcium_tau * std::exp((-calcium_tau * std::log(calcium_tau / decay_tau)) / (calcium_tau - decay_tau)));
  }; // Checked in python, this is valid
  // calciumExtrusion is a rate, not a tau
  this->constants.preCalciumFluxFactor  = amplitudeNormalization(prespikeCalcium, 1 / this->constants.calciumExtrusionCtt, preCalciumDecayTau);
  this->constants.postCalciumFluxFactor = amplitudeNormalization(postspikeCalcium, 1 / this->constants.calciumExtrusionCtt, postCalciumDecayTau);
}

void MACRbPModel::CheckParameters(const std::vector<FileEntry> &parameters) {
  BranchedMorphology::CheckParameters(parameters);
  for (auto &[parameterName, parameterValues] : parameters) {
    if (parameterName.find("kinasesTotal") != std::string::npos && (this->constants.kinasesTotal != std::stod(parameterValues.at(0)))) {
      throw "kinasesTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("phosphatasesTotal") != std::string::npos && this->constants.calcineurinTotal != std::stod(parameterValues.at(0))) {
      throw "phosphatasesTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("calmodulinTotal") != std::string::npos && this->constants.calmodulinTotal != std::stod(parameterValues.at(0))) {
      throw "calmodulinTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("neurograninTotal") != std::string::npos && this->constants.neurograninTotal != std::stod(parameterValues.at(0))) {
      throw "neurograninTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kOne") != std::string::npos &&
               this->constants.reaction1Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kOne was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kTwo") != std::string::npos &&
               this->constants.reaction2Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kTwo was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kThree") != std::string::npos &&
               this->constants.reaction3Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kThree was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kFour") != std::string::npos &&
               this->constants.reaction4Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kFour was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kFive") != std::string::npos &&
               this->constants.reaction5Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kFive was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kSix") != std::string::npos &&
               this->constants.reaction6Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kSix was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kSeven") != std::string::npos &&
               this->constants.reaction7Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kSeven was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kEight") != std::string::npos &&
               this->constants.reaction8Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kEight was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kNine") != std::string::npos &&
               this->constants.reaction9Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kNine was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kTen") != std::string::npos &&
               this->constants.reaction10Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kTen was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kEleven") != std::string::npos &&
               this->constants.reaction11Ctt / constants.resourceConversionFct != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kEleven was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kTwelve") != std::string::npos &&
               this->constants.reaction12Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kTwelve was not consistent in plasticity model parameters.";
    } else if (parameterName.find("resourceDiffusion") != std::string::npos &&
               this->constants.resourceDiffusionFct != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2)) {
      throw "resourceDiffusion was not consistent in plasticity model parameters.";
    } else if (parameterName.find("calciumDiffusion") != std::string::npos &&
               this->constants.caDiffusionFct != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2)) {
      throw "calciumDiffusion was not consistent in plasticity model parameters.";
    } else if (parameterName.find("calciumExtrusion") != std::string::npos &&
               this->constants.calciumExtrusionCtt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "calciumExtrusion was not consistent in plasticity model parameters.";
    } else if (parameterName.find("initialWeight") != std::string::npos && this->constants.initialWeight != std::stod(parameterValues.at(0))) {
      throw "initialWeight was not consistent in plasticity model parameters.";
    } else if (parameterName.find("availResourcesRatio") != std::string::npos && this->availResourcesRatio != std::stod(parameterValues.at(0))) {
      throw "availResourcesRatio was not consistent in plasticity model parameters.";
    } else if (parameterName.find("prespikeCalcium") != std::string::npos && this->prespikeCalcium != std::stod(parameterValues.at(0))) {
      throw "prespikeCalcium was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postspikeCalcium") != std::string::npos && this->postspikeCalcium != std::stod(parameterValues.at(0))) {
      throw "postspikeCalcium was not consistent in plasticity model parameters.";
    } else if (parameterName.find("preCalciumRiseTau") != std::string::npos &&
               this->preCalciumRiseTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("preCalciumDecayTau") != std::string::npos &&
               this->preCalciumDecayTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postCalciumRiseTau") != std::string::npos &&
               this->postCalciumRiseTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postCalciumDecayTau") != std::string::npos &&
               this->postCalciumDecayTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("nonlinearNMDAFactor") != std::string::npos &&
               this->constants.nonlinearFactorNMDA != (std::stod(parameterValues.at(0)))) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("resourceConversionFct") != std::string::npos &&
               this->constants.resourceConversionFct != (std::stod(parameterValues.at(0)))) {
      throw "#mV*uM^-1 Conversion factor from nM to mV";

    } else if (parameterName.find("calciumBasal") != std::string::npos && this->calciumBasal != (std::stod(parameterValues.at(0)))) {
      throw "calciumBasal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("newtonIterations") != std::string::npos &&
               this->constants.newtonIterations != (std::stod(parameterValues.at(0)))) {
      throw "newtonIterations iterations was not consistent in plasticity model parameters.";
    }
  }
}

void MACRbPModel::SaveParameters(std::ofstream &wParameterFile, std::string neuronIdentificator) const {
  BranchedMorphology::SaveParameters(wParameterFile, neuronIdentificator);
  // Parameter outputs should not be written in std::to_string unless you are confident that the value will not be too small for the I/O to work as
  // intended, or you will lose information in the Parameter.txt file
  wParameterFile << neuronIdentificator << "kinasesTotal\t\t" << (this->constants.kinasesTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of CaMKII in the synapse spine\n";
  wParameterFile << neuronIdentificator << "phosphatasesTotal\t\t" << (this->constants.calcineurinTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of calcineurin in the synapse spine\n";
  wParameterFile << neuronIdentificator << "calmodulinTotal\t\t" << (this->constants.calmodulinTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of CaM in the synapse spine\n";
  wParameterFile << neuronIdentificator << "neurograninTotal\t\t" << (this->constants.neurograninTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of neurogranin in the synapse spine\n";

  wParameterFile << neuronIdentificator << "newtonIterations\t\t" << this->constants.newtonIterations << " #_";
  wParameterFile << "\t"
                 << "#Newton iterations for the neurogranin equilibrium\n";
  wParameterFile << neuronIdentificator << "kOne\t\t" << (this->constants.reaction1Ctt / infoGlobal->dtTimestep) << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From CaM and Ng to CaMNg\n";
  wParameterFile << neuronIdentificator << "kTwo\t\t" << (this->constants.reaction2Ctt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#From CaMNg to CaM and Ng\n";
  wParameterFile << neuronIdentificator << "kThree\t\t" << (this->constants.reaction3Ctt / infoGlobal->dtTimestep) << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From CaM and Ca to active CaM\n";
  wParameterFile << neuronIdentificator << "kFour\t\t" << (this->constants.reaction4Ctt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#From active CaM to CaM and Ca\n";
  wParameterFile << neuronIdentificator << "kFive\t\t" << (this->constants.reaction5Ctt / infoGlobal->dtTimestep) << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From N inactive to N active, using CaM\n";
  wParameterFile << neuronIdentificator << "kSix\t\t" << (this->constants.reaction6Ctt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#From N active to N inactive, using CaM\n";
  wParameterFile << neuronIdentificator << "kSeven\t\t" << (this->constants.reaction7Ctt / infoGlobal->dtTimestep) << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From K bound to CaM to phosphorylated K\n";
  wParameterFile << neuronIdentificator << "kEight\t\t" << (this->constants.reaction8Ctt / infoGlobal->dtTimestep) << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From phosphorylated K to K bound to CaM\n";
  wParameterFile << neuronIdentificator << "kNine\t\t" << (this->constants.reaction9Ctt / infoGlobal->dtTimestep) << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From K inactive to K active, using CaM\n";
  wParameterFile << neuronIdentificator << "kTen\t\t" << (this->constants.reaction10Ctt / infoGlobal->dtTimestep) << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From K active to K inactive, using CaM\n";
  wParameterFile << neuronIdentificator << "kEleven\t\t" << (this->constants.reaction11Ctt / constants.resourceConversionFct / infoGlobal->dtTimestep)
                 << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From resources to weight\n";
  wParameterFile << neuronIdentificator << "kTwelve\t\t" << (this->constants.reaction12Ctt / infoGlobal->dtTimestep) << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From weight to resources\n";

  wParameterFile << neuronIdentificator << "resourceDiffusion\t"
                 << (this->constants.resourceDiffusionFct * std::pow(synapticGap, 2) / infoGlobal->dtTimestep) << " #um^{2} s^{-1}";
  wParameterFile << "\t"
                 << "#Diffusion constant of resources used to increase synapse size\n";
  wParameterFile << neuronIdentificator << "calciumDiffusion\t"
                 << (this->constants.caDiffusionFct * std::pow(synapticGap, 2) / infoGlobal->dtTimestep) << " #um^{2} s^{-1}";
  wParameterFile << "\t"
                 << "#Diffusion constant of free calcium in the branch\n";

  wParameterFile << neuronIdentificator << "initialWeight\t\t" << (this->constants.initialWeight) << " #mV/spike";
  wParameterFile << "\t"
                 << "#Initial weight value\n";
  wParameterFile << neuronIdentificator << "availResourcesRatio\t\t" << (this->availResourcesRatio) << "#-";
  wParameterFile << "\t"
                 << "#Initial amount of resources for each segment of the dendrite\n";
  wParameterFile << neuronIdentificator << "prespikeCalcium\t\t" << (this->prespikeCalcium) << " #nM/spike";
  wParameterFile << "\t"
                 << "#Total influx of calcium with a prespike\n";
  wParameterFile << neuronIdentificator << "postspikeCalcium\t\t" << (this->postspikeCalcium) << " #nM/spike";
  wParameterFile << "\t"
                 << "#Total influx of calcium with a postspike\n";
  wParameterFile << neuronIdentificator << "preCalciumRiseTau\t\t" << (static_cast<double>(this->preCalciumRiseTau) * infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "preCalciumDecayTau\t\t" << (static_cast<double>(this->preCalciumDecayTau) * infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "postCalciumRiseTau\t\t" << (static_cast<double>(this->postCalciumRiseTau) * infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "postCalciumDecayTau\t\t" << (static_cast<double>(this->postCalciumDecayTau) * infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "nonlinearNMDAFactor\t\t" << (static_cast<double>(this->constants.nonlinearFactorNMDA)) << " #-";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";

  wParameterFile << neuronIdentificator << "calciumBasal\t\t" << (calciumBasal) << " #uM";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "calciumExtrusion\t\t" << (this->constants.calciumExtrusionCtt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#Extrusion rate for free calcium\n";
  wParameterFile << neuronIdentificator << "resourceConversionFct\t\t" << (constants.resourceConversionFct) << " #mV spike^{-1} uM^{-1}";
  wParameterFile << "\t"
                 << "#Conversion from resource molar units (uM) to dmV/spike\n";
}
// Remove from here
MACRbPSynapseSpine MACRbPModel::ComputeSteadyState(MACRbPSynapseSpine &spine) {
  // Set values
  spine.calciumFree        = calciumBasal;
  spine.resourcesAvailable = constants.initialResources;
  spine.weight             = constants.initialWeight;
  // Time stamp
  auto steady_state_begin = std::chrono::high_resolution_clock::now();

  // Run loop
  while (std::abs(spine.resourcesAvailable - spine.resourcesOldStep) > std::numeric_limits<double>::epsilon()) {
    // Calcium diffusion
    spine.calciumFree += constants.calciumInfluxBasal - spine.calciumFree * constants.calciumExtrusionCtt;
    spine.resourcesOldStep = spine.resourcesAvailable;

    double freeNg{constants.neurograninTotal - spine.calmodulinNeurogranin}, epsilon{};
    for (int i = 0; i < constants.newtonIterations; i++) {
      epsilon = (constants.reaction1234Ctt * freeNg * spine.calmodulinActive - spine.calmodulinNeurogranin * std::pow(spine.calciumFree, 2)) /
                (constants.reaction1234Ctt * (freeNg + spine.calmodulinActive) + 4 * spine.calmodulinNeurogranin * spine.calciumFree +
                 std::pow(spine.calciumFree, 2));
      freeNg -= epsilon;
      spine.calmodulinActive -= epsilon;
      spine.calmodulinNeurogranin += epsilon;
      spine.calciumFree += 2 * epsilon;
    }

    double kDot{}, nDot{}, kPDot{}, wDot{};
    nDot = constants.reaction5Ctt * (constants.calcineurinTotal - spine.calcineurinActive) * spine.calmodulinActive -
           constants.reaction6Ctt * spine.calcineurinActive;
    spine.calcineurinActive += nDot;
    spine.calmodulinActive -= nDot;
    // 5th kinase autophosphorylation
    kPDot = constants.reaction7Ctt * spine.kinasesCaM * (spine.kinasesCaM + spine.kinasesPhospho) -
            constants.reaction8Ctt * spine.calcineurinActive * spine.kinasesPhospho;
    spine.kinasesPhospho += kPDot;
    spine.kinasesCaM -= kPDot;
    // 6th kinase activation via CaM -kP
    kDot = constants.reaction9Ctt * (constants.kinasesTotal - spine.kinasesCaM - spine.kinasesPhospho) * spine.calmodulinActive -
           constants.reaction10Ctt * spine.kinasesCaM;
    spine.kinasesCaM += kDot;
    spine.calmodulinActive -= kDot;
    // ORDERING
    //  7th active CaM consumption by Ndot and Kdot (not kPdot)
    // spine.calmodulinActive -= nDot + kDot;//We are not doing this because this could mean negative active calmodulins
    // 8th Change in synapse spine size/weight
    wDot = constants.reaction11Ctt * spine.resourcesAvailable * (spine.kinasesCaM + spine.kinasesPhospho) -
           constants.reaction12Ctt * spine.weight * spine.calcineurinActive;
    spine.weight += wDot;
    // 9th Consumption of resources by weight change
    spine.resourcesAvailable -=
        (wDot) / constants.resourceConversionFct; // This should be the case for converting the dmV/spike to concentration of resources
    // Store equilibrium in a file? Together with time spent in the calculation? Output only to terminal the time?
  }

  auto steady_state_end = std::chrono::high_resolution_clock::now();
  auto ss_time          = std::chrono::duration_cast<std::chrono::seconds>(steady_state_end - steady_state_begin);
  std::cout << "Time to reach steady state: " << ss_time.count() << '\n';
  return spine;
}

int MACRbPModel::CreateBranch(std::vector<int> anteriorBranches) {
  // After the cheap fix, we do what we used to do
  int branchId{this->GenerateBranchId()};
  if (!anteriorBranches.empty()) {
    this->MACRbPBranches.push_back(MACRbPBranch(anteriorBranches, this->synapticGap, this->branchLength, branchId, constants,
                                                steady_state_spine)); // This vector should be sorted by ID by default (tested).
    this->branches.push_back(static_cast<Branch *>(&this->MACRbPBranches.back()));
  } else {
    int branchId{this->GenerateBranchId()};
    this->MACRbPBranches.push_back(MACRbPBranch(this->synapticGap, this->branchLength, branchId, constants,
                                                steady_state_spine)); // This vector should be sorted by ID by default (tested).
    this->branches.push_back(static_cast<Branch *>(&this->MACRbPBranches.back()));
  }
  return branchId;
}

void MACRbPModel::Advect() {
  // if (!postSpiked){
  //     TStepModded=infoGlobal->timeStep%preSpikeDelaySteps;
  // }
  // TStepModded++;//Here we could assume that the timestep goes hand in hand with this integer in terms of increments.
  // TStepModded%=preSpikeDelaySteps;
  // TStepInput=TStepModded+preSpikeDelaySteps+1;
  // TStepInput%=preSpikeDelaySteps;
  for (MACRbPBranch &branch : MACRbPBranches) {
    branch.Advect();
    // branch.Advect(TStepModded);
  }
  // TStepModded++;//This is for the postspike calcium influx next timestep
}

void MACRbPModel::RecordPostSpike() {
  // TStepModded=infoGlobal->timeStep%preSpikeDelaySteps;
  this->totalPostSpikes++;
  // this->postSpiked = true;
  for (MACRbPBranch &branch : MACRbPBranches) {
    branch.PostSpikeCalciumFlux();
  }
  // for (CaDiffusionBranch& branch: caDiffBranches){
  // 	std::transform(PAR_UNSEQ,branch.waitingMatrix.at(TStepModded).begin(), branch.waitingMatrix.at(TStepModded).end(),
  // branch.waitingMatrix.at(TStepModded).begin(),std::bind(std::plus<double>(),std::placeholders::_1, postspikeCalcium));
  // }
}

void MACRbPModel::RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) {
  std::lock_guard<std::mutex> _locked(_presynSpikeMutex);
  // CaResSynapseSpine& synapseSpine = *caResSpines.at(spikedSpineId);
  // CaDiffusionBranch&  branch       = caDiffBranches.at(synapseSpine.branchId);
  // int  branchSpinePosition{synapseSpine.branchPositionId};
  static_cast<MACRbpSpinePtr>(spinePtr)->preTransientIncrease++;
  // branch.waitingMatrix.at(TStepInput).at(branchSpinePosition)+=prespikeCalcium;
  this->totalPreSpikes++;
}

void MACRbPModel::PostConnectSetUp() {
  for (MACRbPBranch &branch : MACRbPBranches) {
    branch.PostConnectSetUp(branchedSpineData);
  }
}

BaseSpinePtr MACRbPModel::AllocateNewSynapse(BranchTargeting &branchTarget) {
  int branch{AllocateBranch(branchTarget)};
  int position{PopSynapseSlotFromBranch(branchTarget)};
  // caDiffBranches.at(branch).CaDiffSpines.push_back(CaResSynapseSpine(kinasesTotal, calcineurinTotal, constants.initialWeight));
  MACRbpSpinePtr newSpine = &MACRbPBranches.at(branch).MACRbPspines.at(position);

  // this->weightsSum += newSynapse->GetWeight();
  newSpine->idInMorpho = (this->baseSpineData.size()); // this->spineIdGenerator++
  newSpine->weight     = constants.initialWeight;
  // Branch
  newSpine->branchPositionId = (position);
  newSpine->branchId         = (branch);

  newSpine->connected = true;

  branches.at(branch)->synapseSlotClosedIndex.push_back(position); // Do we really need this?

  // Storage (other)
  this->baseSpineData.push_back(static_cast<BaseSpinePtr>(newSpine));
  this->branchedSpineData.push_back(static_cast<BranchedSpinePtr>(newSpine));
  // this->MACRbPSpines.push_back(newSpine);

  return this->baseSpineData.back();
}

std::vector<double> MACRbPModel::GetOverallSynapticProfile() const {
  // What could we put here? Weight disparity(deviation of weight between neighbouring spines)? Avg weight might help
  std::vector<double> dataArray(4);
  size_t              noSpines{this->baseSpineData.size()};
  double              total_weight{};
  for (const MACRbPBranch &branch : MACRbPBranches) {
    total_weight += branch.GetTotalWeight();
  }

  dataArray.at(0) = total_weight / noSpines;
  dataArray.at(1) = this->totalPostSpikes;
  dataArray.at(2) = this->totalPreSpikes;
  dataArray.at(3) = noSpines;
  return dataArray;
}

std::string MACRbPModel::GetOverallSynapticProfileHeaderInfo() const {
  return std::string("{<average weight>, <total post spikes> ,<total pre spikes>, <total number of spines>}");
}

std::vector<double> MACRbPModel::GetSteadyStateData() const {
  return std::vector<double>();
}

std::vector<std::string> MACRbPModel::GetSteadyStateVarNames() const {
  return std::vector<std::string>{
      "",
  };
}
