#include "./MACRbPModel.hpp"

MACRbPModel::MACRbPModel(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters, MACRbPSynapseSpine &spine):
    BranchedMorphology(infoGlobal, morphologyParameters) {
  if (spine.calmodulinActive == 0) {
    ComputeSteadyState(spine);
  }
  steady_state_spine = spine;
  // constant_parameters.newtonIterations = 2;
}

void MACRbPModel::LoadParameters(const std::vector<FileEntry> &morphologyParameters) {
  BranchedMorphology::LoadParameters(morphologyParameters);
  for (auto &[parameterName, parameterValues] : morphologyParameters) {
    if (parameterName.find("kinasesTotal") != std::string::npos) {
      this->runtime_constants.kinasesTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("phosphatasesTotal") != std::string::npos) {
      this->runtime_constants.calcineurinTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("calmodulinTotal") != std::string::npos) {
      this->runtime_constants.calmodulinTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("neurograninTotal") != std::string::npos) {
      this->runtime_constants.neurograninTotal = std::stod(parameterValues.at(0));
    } else if (parameterName.find("kOne") != std::string::npos) {
      this->constant_parameters.reaction1Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kTwo") != std::string::npos) {
      this->constant_parameters.reaction2Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kThree") != std::string::npos) {
      this->constant_parameters.reaction3Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kFour") != std::string::npos) {
      this->constant_parameters.reaction4Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kFive") != std::string::npos) {
      this->runtime_constants.reaction5Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kSix") != std::string::npos) {
      this->runtime_constants.reaction6Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kSeven") != std::string::npos) {
      this->runtime_constants.reaction7Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kEight") != std::string::npos) {
      this->runtime_constants.reaction8Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kNine") != std::string::npos) {
      this->runtime_constants.reaction9Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kTen") != std::string::npos) {
      this->runtime_constants.reaction10Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kEleven") != std::string::npos) {
      this->runtime_constants.reaction11Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("kTwelve") != std::string::npos) {
      this->runtime_constants.reaction12Ctt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("resourceDiffusion") != std::string::npos) {
      this->runtime_constants.resourceDiffusionFct = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2);
    } else if (parameterName.find("calciumDiffusion") != std::string::npos) {
      this->runtime_constants.caDiffusionFct = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2);
    } else if (parameterName.find("calciumExtrusion") != std::string::npos) {
      this->runtime_constants.calciumExtrusionCtt = std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep;
    } else if (parameterName.find("initialWeight") != std::string::npos) {
      this->constant_parameters.initialWeight = std::stod(parameterValues.at(0));
    } else if (parameterName.find("initialResources") != std::string::npos) {
      this->constant_parameters.availResourcesRatio = std::stod(parameterValues.at(0));
    } else if (parameterName.find("prespikeCalcium") != std::string::npos) {
      this->constant_parameters.prespikeCalcium = std::stod(parameterValues.at(0));
    } else if (parameterName.find("postspikeCalcium") != std::string::npos) {
      this->constant_parameters.postspikeCalcium = std::stod(parameterValues.at(0));
      // } else if (parameterName.find("preCalciumRiseTau") != std::string::npos) {
      //   this->preCalciumRiseTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("preCalciumDecayTau") != std::string::npos) {
      this->constant_parameters.preCalciumDecayTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
      // } else if (parameterName.find("postCalciumRiseTau") != std::string::npos) {
      //   this->postCalciumRiseTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("postCalciumDecayTau") != std::string::npos) {
      this->constant_parameters.postCalciumDecayTau = (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep);
    } else if (parameterName.find("nonlinearNMDAFactor") != std::string::npos) {
      this->runtime_constants.nonlinearFactorNMDA = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("calciumBasal") != std::string::npos) {
      this->constant_parameters.calciumBasal = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("resourceConversionFct") != std::string::npos) {
      this->constant_parameters.resourceConversionFct = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("newtonIterations") != std::string::npos) {
      this->runtime_constants.newtonIterations = (std::stoi(parameterValues.at(0)));
    }
  }
  PreCalculateParameters();
}

void MACRbPModel::PreCalculateParameters() {
  this->runtime_constants.reaction1234Ctt = this->constant_parameters.reaction1Ctt * this->constant_parameters.reaction4Ctt /
                                            (this->constant_parameters.reaction2Ctt * this->constant_parameters.reaction3Ctt);
  this->runtime_constants.calciumInfluxBasal = this->constant_parameters.calciumBasal * runtime_constants.calciumExtrusionCtt;
  this->constant_parameters.initialResources = this->constant_parameters.initialWeight / this->constant_parameters.resourceConversionFct /
                                               this->constant_parameters.availResourcesRatio;  // REVIEW
  // this->constant_parameters.reaction11Ctt *= constant_parameters.resourceConversionFct; // REVIEW
  // this->constant_parameters.preCalciumRiseRate   = 1 / preCalciumRiseTau;
  this->runtime_constants.preCalciumDecayRate = 1 / this->constant_parameters.preCalciumDecayTau;
  // this->constant_parameters.postCalciumRiseRate  = 1 / postCalciumRiseTau;
  this->runtime_constants.postCalciumDecayRate = 1 / this->constant_parameters.postCalciumDecayTau;

  auto amplitudeNormalization = [](double amplitude, double calcium_tau,
                                   double decay_tau) {  // The taus need to be in the proper time unit
    return amplitude / (calcium_tau * std::exp((-calcium_tau * std::log(calcium_tau / decay_tau)) / (calcium_tau - decay_tau)));
  };  // Checked in python, this is valid
  // calciumExtrusion is a rate, not a tau
  this->runtime_constants.preCalciumFluxFactor =
    amplitudeNormalization(this->constant_parameters.prespikeCalcium, 1 / this->runtime_constants.calciumExtrusionCtt,
                           this->constant_parameters.preCalciumDecayTau);
  this->runtime_constants.postCalciumFluxFactor =
    amplitudeNormalization(this->constant_parameters.postspikeCalcium, 1 / this->runtime_constants.calciumExtrusionCtt,
                           this->constant_parameters.postCalciumDecayTau);
}

void MACRbPModel::CheckParameters(const std::vector<FileEntry> &parameters) {
  BranchedMorphology::CheckParameters(parameters);
  for (auto &[parameterName, parameterValues] : parameters) {
    if (parameterName.find("kinasesTotal") != std::string::npos &&
        (this->runtime_constants.kinasesTotal != std::stod(parameterValues.at(0)))) {
      throw "kinasesTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("phosphatasesTotal") != std::string::npos &&
               this->runtime_constants.calcineurinTotal != std::stod(parameterValues.at(0))) {
      throw "phosphatasesTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("calmodulinTotal") != std::string::npos &&
               this->runtime_constants.calmodulinTotal != std::stod(parameterValues.at(0))) {
      throw "calmodulinTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("neurograninTotal") != std::string::npos &&
               this->runtime_constants.neurograninTotal != std::stod(parameterValues.at(0))) {
      throw "neurograninTotal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kOne") != std::string::npos &&
               this->constant_parameters.reaction1Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kOne was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kTwo") != std::string::npos &&
               this->constant_parameters.reaction2Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kTwo was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kThree") != std::string::npos &&
               this->constant_parameters.reaction3Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kThree was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kFour") != std::string::npos &&
               this->constant_parameters.reaction4Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kFour was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kFive") != std::string::npos &&
               this->runtime_constants.reaction5Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kFive was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kSix") != std::string::npos &&
               this->runtime_constants.reaction6Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kSix was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kSeven") != std::string::npos &&
               this->runtime_constants.reaction7Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kSeven was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kEight") != std::string::npos &&
               this->runtime_constants.reaction8Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kEight was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kNine") != std::string::npos &&
               this->runtime_constants.reaction9Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kNine was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kTen") != std::string::npos &&
               this->runtime_constants.reaction10Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kTen was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kEleven") != std::string::npos &&
               this->runtime_constants.reaction11Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kEleven was not consistent in plasticity model parameters.";
    } else if (parameterName.find("kTwelve") != std::string::npos &&
               this->runtime_constants.reaction12Ctt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "kTwelve was not consistent in plasticity model parameters.";
    } else if (parameterName.find("resourceDiffusion") != std::string::npos &&
               this->runtime_constants.resourceDiffusionFct !=
                 std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2)) {
      throw "resourceDiffusion was not consistent in plasticity model parameters.";
    } else if (parameterName.find("calciumDiffusion") != std::string::npos &&
               this->runtime_constants.caDiffusionFct !=
                 std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep / std::pow(synapticGap, 2)) {
      throw "calciumDiffusion was not consistent in plasticity model parameters.";
    } else if (parameterName.find("calciumExtrusion") != std::string::npos &&
               this->runtime_constants.calciumExtrusionCtt != std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep) {
      throw "calciumExtrusion was not consistent in plasticity model parameters.";
    } else if (parameterName.find("initialWeight") != std::string::npos &&
               this->constant_parameters.initialWeight != std::stod(parameterValues.at(0))) {
      throw "initialWeight was not consistent in plasticity model parameters.";
    } else if (parameterName.find("availResourcesRatio") != std::string::npos &&
               this->constant_parameters.availResourcesRatio != std::stod(parameterValues.at(0))) {
      throw "availResourcesRatio was not consistent in plasticity model parameters.";
    } else if (parameterName.find("prespikeCalcium") != std::string::npos &&
               this->constant_parameters.prespikeCalcium != std::stod(parameterValues.at(0))) {
      throw "prespikeCalcium was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postspikeCalcium") != std::string::npos &&
               this->constant_parameters.postspikeCalcium != std::stod(parameterValues.at(0))) {
      throw "postspikeCalcium was not consistent in plasticity model parameters.";
      // } else if (parameterName.find("preCalciumRiseTau") != std::string::npos &&
      //            this->preCalciumRiseTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      //   throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("preCalciumDecayTau") != std::string::npos &&
               this->constant_parameters.preCalciumDecayTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
      // } else if (parameterName.find("postCalciumRiseTau") != std::string::npos &&
      //            this->postCalciumRiseTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      //   throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postCalciumDecayTau") != std::string::npos &&
               this->constant_parameters.postCalciumDecayTau != (std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("nonlinearNMDAFactor") != std::string::npos &&
               this->runtime_constants.nonlinearFactorNMDA != (std::stod(parameterValues.at(0)))) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("resourceConversionFct") != std::string::npos &&
               this->constant_parameters.resourceConversionFct != (std::stod(parameterValues.at(0)))) {
      throw "#mV*uM^-1 Conversion factor from nM to mV";

    } else if (parameterName.find("calciumBasal") != std::string::npos &&
               this->constant_parameters.calciumBasal != (std::stod(parameterValues.at(0)))) {
      throw "calciumBasal was not consistent in plasticity model parameters.";
    } else if (parameterName.find("newtonIterations") != std::string::npos &&
               this->runtime_constants.newtonIterations != (std::stod(parameterValues.at(0)))) {
      throw "newtonIterations iterations was not consistent in plasticity model parameters.";
    }
  }
}

void MACRbPModel::SaveParameters(std::ofstream &wParameterFile, std::string neuronIdentificator) const {
  BranchedMorphology::SaveParameters(wParameterFile, neuronIdentificator);
  // Parameter outputs should not be written in std::to_string unless you are confident that the value will not be too
  // small for the I/O to work as intended, or you will lose information in the Parameter.txt file
  wParameterFile << neuronIdentificator << "kinasesTotal\t\t" << (this->runtime_constants.kinasesTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of CaMKII in the synapse spine\n";
  wParameterFile << neuronIdentificator << "phosphatasesTotal\t\t" << (this->runtime_constants.calcineurinTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of calcineurin in the synapse spine\n";
  wParameterFile << neuronIdentificator << "calmodulinTotal\t\t" << (this->runtime_constants.calmodulinTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of CaM in the synapse spine\n";
  wParameterFile << neuronIdentificator << "neurograninTotal\t\t" << (this->runtime_constants.neurograninTotal) << " #uM/spine";
  wParameterFile << "\t"
                 << "#Concentration of neurogranin in the synapse spine\n";

  wParameterFile << neuronIdentificator << "newtonIterations\t\t" << this->runtime_constants.newtonIterations << " #_";
  wParameterFile << "\t"
                 << "#Newton iterations for the neurogranin equiconstant_parameters";
  wParameterFile << neuronIdentificator << "kOne\t\t" << (this->constant_parameters.reaction1Ctt / infoGlobal->dtTimestep)
                 << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From CaM and Ng to CaMNg\n";
  wParameterFile << neuronIdentificator << "kTwo\t\t" << (this->constant_parameters.reaction2Ctt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#From CaMNg to CaM and Ng\n";
  wParameterFile << neuronIdentificator << "kThree\t\t" << (this->constant_parameters.reaction3Ctt / infoGlobal->dtTimestep)
                 << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From CaM and Ca to active CaM\n";
  wParameterFile << neuronIdentificator << "kFour\t\t" << (this->constant_parameters.reaction4Ctt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#From active CaM to CaM and Ca\n";
  wParameterFile << neuronIdentificator << "kFive\t\t" << (this->runtime_constants.reaction5Ctt / infoGlobal->dtTimestep)
                 << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From N inactive to N active, using CaM\n";
  wParameterFile << neuronIdentificator << "kSix\t\t" << (this->runtime_constants.reaction6Ctt / infoGlobal->dtTimestep) << " #s^{-1}";
  wParameterFile << "\t"
                 << "#From N active to N inactive, using CaM\n";
  wParameterFile << neuronIdentificator << "kSeven\t\t" << (this->runtime_constants.reaction7Ctt / infoGlobal->dtTimestep)
                 << " #uM^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From K bound to CaM to phosphorylated K\n";
  wParameterFile << neuronIdentificator << "kEight\t\t" << (this->runtime_constants.reaction8Ctt / infoGlobal->dtTimestep)
                 << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From phosphorylated K to K bound to CaM\n";
  wParameterFile << neuronIdentificator << "kNine\t\t" << (this->runtime_constants.reaction9Ctt / infoGlobal->dtTimestep)
                 << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From K inactive to K active, using CaM\n";
  wParameterFile << neuronIdentificator << "kTen\t\t" << (this->runtime_constants.reaction10Ctt / infoGlobal->dtTimestep)
                 << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From K active to K inactive, using CaM\n";
  wParameterFile << neuronIdentificator << "kEleven\t\t" << ((this->runtime_constants.reaction11Ctt) / infoGlobal->dtTimestep)
                 << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From resources to weight\n";
  wParameterFile << neuronIdentificator << "kTwelve\t\t" << (this->runtime_constants.reaction12Ctt / infoGlobal->dtTimestep)
                 << " #u M^{-1} s^{-1}";
  wParameterFile << "\t"
                 << "#From weight to resources\n";

  wParameterFile << neuronIdentificator << "resourceDiffusion\t"
                 << (this->runtime_constants.resourceDiffusionFct * std::pow(synapticGap, 2) / infoGlobal->dtTimestep) << " #um^{2} s^{-1}";
  wParameterFile << "\t"
                 << "#Diffusion constant of resources used to increase synapse size\n";
  wParameterFile << neuronIdentificator << "calciumDiffusion\t"
                 << (this->runtime_constants.caDiffusionFct * std::pow(synapticGap, 2) / infoGlobal->dtTimestep) << " #um^{2} s^{-1}";
  wParameterFile << "\t"
                 << "#Diffusion constant of free calcium in the branch\n";

  wParameterFile << neuronIdentificator << "initialWeight\t\t" << (this->constant_parameters.initialWeight) << " #mV/spike";
  wParameterFile << "\t"
                 << "#Initial weight value\n";
  wParameterFile << neuronIdentificator << "availResourcesRatio\t\t" << (this->constant_parameters.availResourcesRatio) << "#-";
  wParameterFile << "\t"
                 << "#Initial amount of resources for each segment of the dendrite\n";
  wParameterFile << neuronIdentificator << "prespikeCalcium\t\t" << (this->constant_parameters.prespikeCalcium) << " #nM/spike";
  wParameterFile << "\t"
                 << "#Total influx of calcium with a prespike\n";
  wParameterFile << neuronIdentificator << "postspikeCalcium\t\t" << (this->constant_parameters.postspikeCalcium) << " #nM/spike";
  wParameterFile << "\t"
                 << "#Total influx of calcium with a postspike\n";
  // wParameterFile << neuronIdentificator << "preCalciumRiseTau\t\t" << (static_cast<double>(this->preCalciumRiseTau) *
  // infoGlobal->dtTimestep)
  //                << " #secs.";
  // wParameterFile << "\t"
  //                << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "preCalciumDecayTau\t\t"
                 << (static_cast<double>(this->constant_parameters.preCalciumDecayTau) * infoGlobal->dtTimestep) << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  // wParameterFile << neuronIdentificator << "postCalciumRiseTau\t\t" << (static_cast<double>(this->postCalciumRiseTau)
  // * infoGlobal->dtTimestep)
  //                << " #secs.";
  // wParameterFile << "\t"
  //                << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "postCalciumDecayTau\t\t"
                 << (static_cast<double>(this->constant_parameters.postCalciumDecayTau) * infoGlobal->dtTimestep) << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "nonlinearNMDAFactor\t\t" << (static_cast<double>(this->runtime_constants.nonlinearFactorNMDA))
                 << " #-";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";

  wParameterFile << neuronIdentificator << "calciumBasal\t\t" << (constant_parameters.calciumBasal) << " #uM";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "calciumExtrusion\t\t" << (this->runtime_constants.calciumExtrusionCtt / infoGlobal->dtTimestep)
                 << " #s^{-1}";
  wParameterFile << "\t"
                 << "#Extrusion rate for free calcium\n";
  wParameterFile << neuronIdentificator << "resourceConversionFct\t\t" << (constant_parameters.resourceConversionFct)
                 << " #mV spike^{-1} uM^{-1}";
  wParameterFile << "\t"
                 << "#Conversion from resource molar units (uM) to dmV/spike\n";
}
// Remove from here
MACRbPSynapseSpine MACRbPModel::ComputeSteadyState(MACRbPSynapseSpine &spine) {
  // Set values
  spine.calciumFree        = constant_parameters.calciumBasal;
  spine.resourcesAvailable = constant_parameters.initialResources;
  spine.weight = constant_parameters.initialWeight / constant_parameters.resourceConversionFct;  // Weight will be treated as uMs
  spine.couplingStrength /= constant_parameters.resourceConversionFct;
  // Time stamp
  auto steady_state_begin = std::chrono::high_resolution_clock::now();

  // Run loop
  while (std::abs(spine.resourcesAvailable - spine.resourcesOldStep) > std::numeric_limits<double>::epsilon()) {
    // Calcium diffusion
    spine.calciumFree += runtime_constants.calciumInfluxBasal - spine.calciumFree * runtime_constants.calciumExtrusionCtt;
    spine.resourcesOldStep = spine.resourcesAvailable;

    double freeNg { runtime_constants.neurograninTotal - spine.calmodulinNeurogranin }, epsilon {};
    for (int i = 0; i < runtime_constants.newtonIterations; i++) {
      epsilon = (runtime_constants.reaction1234Ctt * freeNg * spine.calmodulinActive -
                 spine.calmodulinNeurogranin * std::pow(spine.calciumFree, 2)) /
                (runtime_constants.reaction1234Ctt * (freeNg + spine.calmodulinActive) +
                 4 * spine.calmodulinNeurogranin * spine.calciumFree + std::pow(spine.calciumFree, 2));
      freeNg -= epsilon;
      spine.calmodulinActive -= epsilon;
      spine.calmodulinNeurogranin += epsilon;
      spine.calciumFree += 2 * epsilon;
    }

    double kDot {}, nDot {}, kPDot {}, wDot {};
    nDot = runtime_constants.reaction5Ctt * (runtime_constants.calcineurinTotal - spine.calcineurinActive) * spine.calmodulinActive -
           runtime_constants.reaction6Ctt * spine.calcineurinActive;
    spine.calcineurinActive += nDot;
    spine.calmodulinActive -= nDot;
    // 5th kinase autophosphorylation
    kPDot = runtime_constants.reaction7Ctt * spine.kinasesCaM * (spine.kinasesCaM + spine.kinasesPhospho) -
            runtime_constants.reaction8Ctt * spine.calcineurinActive * spine.kinasesPhospho;
    spine.kinasesPhospho += kPDot;
    spine.kinasesCaM -= kPDot;
    // 6th kinase activation via CaM -kP
    kDot =
      runtime_constants.reaction9Ctt * (runtime_constants.kinasesTotal - spine.kinasesCaM - spine.kinasesPhospho) * spine.calmodulinActive -
      runtime_constants.reaction10Ctt * spine.kinasesCaM;
    spine.kinasesCaM += kDot;
    spine.calmodulinActive -= kDot;
    // ORDERING
    //  7th active CaM consumption by Ndot and Kdot (not kPdot)
    // spine.calmodulinActive -= nDot + kDot;//We are not doing this because this could mean negative active calmodulins
    // 8th Change in synapse spine size/weight
    wDot = runtime_constants.reaction11Ctt * spine.resourcesAvailable * (spine.kinasesCaM + spine.kinasesPhospho) -
           runtime_constants.reaction12Ctt * spine.weight * spine.calcineurinActive;
    // Now wDot is in uM units
    spine.weight += wDot;
    // 9th Consumption of resources by weight change
    spine.resourcesAvailable -= (wDot);
    // / constant_parameters.resourceConversionFct; // This should be the case for converting the dmV/spike to
    // concentration of resources Store equilibrium in a file? Together with time spent in the calculation? Output only
    // to terminal the time?
  }

  auto steady_state_end = std::chrono::high_resolution_clock::now();
  auto ss_time          = std::chrono::duration_cast<std::chrono::seconds>(steady_state_end - steady_state_begin);
  std::cout << "Time to reach steady state: " << ss_time.count() << '\n';
  return spine;
}

int MACRbPModel::CreateBranch(std::vector<int> anteriorBranches) {
  // After the cheap fix, we do what we used to do
  int branchId { this->GenerateBranchId() };
  if (!anteriorBranches.empty()) {
    this->MACRbPBranches.push_back(MACRbPBranch(anteriorBranches, this->synapticGap, this->branchLength, branchId,
                                                steady_state_spine));  // This vector should be sorted by ID by default (tested).
    this->branches.push_back(static_cast<Branch *>(&this->MACRbPBranches.back()));
  } else {
    int branchId { this->GenerateBranchId() };
    this->MACRbPBranches.push_back(MACRbPBranch(this->synapticGap, this->branchLength, branchId,
                                                steady_state_spine));  // This vector should be sorted by ID by default (tested).
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
    branch.Advect(runtime_constants);
    // branch.Advect(TStepModded);
  }
  // TStepModded++;//This is for the postspike calcium influx next timestep
}

void MACRbPModel::RecordPostSpike() {
  // TStepModded=infoGlobal->timeStep%preSpikeDelaySteps;
  this->totalPostSpikes++;
  // this->postSpiked = true;
  for (MACRbPBranch &branch : MACRbPBranches) {
    branch.PostSpikeCalciumFlux(runtime_constants.nonlinearFactorNMDA);
  }
  // for (CaDiffusionBranch& branch: caDiffBranches){
  // 	std::transform(PAR_UNSEQ,branch.waitingMatrix.at(TStepModded).begin(), branch.waitingMatrix.at(TStepModded).end(),
  // branch.waitingMatrix.at(TStepModded).begin(),std::bind(std::plus<double>(),std::placeholders::_1,
  // postspikeCalcium));
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
  int branch { AllocateBranch(branchTarget) };
  int position { PopSynapseSlotFromBranch(branchTarget) };
  // caDiffBranches.at(branch).CaDiffSpines.push_back(CaResSynapseSpine(kinasesTotal, calcineurinTotal,
  // constant_parameters.initialWeight));
  MACRbpSpinePtr newSpine = &MACRbPBranches.at(branch).MACRbPspines.at(position);

  // this->weightsSum += newSynapse->GetWeight();
  newSpine->idInMorpho = (this->baseSpineData.size());  // this->spineIdGenerator++
  // newSpine->weight     = constant_parameters.initialWeight / constant_parameters.resourceConversionFct;//Not
  // necessary, weight is set in the SS computation. This would be an active bug Branch
  newSpine->branchPositionId = (position);
  newSpine->branchId         = (branch);

  newSpine->enabled = true;

  branches.at(branch)->synapseSlotClosedIndex.push_back(position);  // Do we really need this?

  // Storage (other)
  this->baseSpineData.push_back(static_cast<BaseSpinePtr>(newSpine));
  this->branchedSpineData.push_back(static_cast<BranchedSpinePtr>(newSpine));
  // this->MACRbPSpines.push_back(newSpine);

  return this->baseSpineData.back();
}

std::vector<double> MACRbPModel::GetOverallSynapticProfile() const {
  // What could we put here? Weight disparity(deviation of weight between neighbouring spines)? Avg weight might help
  std::vector<double> dataArray(4);
  size_t              noSpines { this->baseSpineData.size() };
  double              total_weight {};
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

std::vector<std::pair<std::string, double>> MACRbPModel::GetSteadyStateData() const {
  return std::vector<std::pair<std::string, double>> {
    { "Calcium ion", steady_state_spine.calciumFree },
    { "Ng-CaM", steady_state_spine.calmodulinNeurogranin },
    { "Free Ng", runtime_constants.neurograninTotal - steady_state_spine.calmodulinNeurogranin },
    { "Ca2CaM", steady_state_spine.calmodulinActive },
    { "Free CaN", runtime_constants.calcineurinTotal - steady_state_spine.calcineurinActive },
    { "Ca2CaM-CaN", runtime_constants.calcineurinTotal - steady_state_spine.calcineurinActive },
    { "Free CaMKII", runtime_constants.kinasesTotal - steady_state_spine.kinasesCaM - steady_state_spine.kinasesPhospho },
    { "Ca2CaM-CaMKII", steady_state_spine.kinasesCaM },
    { "CaMKII-P", steady_state_spine.kinasesPhospho },
    { "AMPAR", steady_state_spine.resourcesAvailable },
    { "Synaptic weight", steady_state_spine.GetWeight() }
  };
}