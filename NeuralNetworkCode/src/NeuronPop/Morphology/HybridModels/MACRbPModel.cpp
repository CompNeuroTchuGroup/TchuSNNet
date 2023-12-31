#include "./MACRbPModel.hpp"

MACRbPModel::MACRbPModel(GlobalSimInfo *infoGlobal) : BranchedMorphology(infoGlobal) {
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
      this->preCalciumRiseTau = (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep);
    } else if (parameterName.find("preCalciumDecayTau") != std::string::npos) {
      this->preCalciumDecayTau = (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep);
    } else if (parameterName.find("postCalciumRiseTau") != std::string::npos) {
      this->postCalciumRiseTau = (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep);
    } else if (parameterName.find("postCalciumDecayTau") != std::string::npos) {
      this->postCalciumDecayTau = (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep);
    } else if (parameterName.find("nonlinearNMDAFactor") != std::string::npos) {
      this->constants.nonlinearFactorNMDA = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("calciumBasal") != std::string::npos) {
      this->calciumBasal = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("resourceConversionFct") != std::string::npos) {
      this->constants.resourceConversionFct = (std::stod(parameterValues.at(0)));
    }
  }
  this->constants.calciumInfluxBasal = calciumBasal * constants.calciumExtrusionCtt;
  this->constants.initialResources   = this->constants.initialWeight / constants.resourceConversionFct / availResourcesRatio;
  this->constants.reaction11Ctt *= constants.resourceConversionFct;
  this->constants.preCalciumRiseRate   = 1 / preCalciumRiseTau;
  this->constants.preCalciumDecayRate  = 1 / preCalciumDecayTau;
  this->constants.postCalciumRiseRate  = 1 / postCalciumRiseTau;
  this->constants.postCalciumDecayRate = 1 / postCalciumDecayTau;

  this->constants.preCalciumFluxFactor =
      prespikeCalcium * std::pow((((1 / preCalciumDecayTau) - (1 / preCalciumRiseTau)) *
                                  (std::pow(preCalciumRiseTau / preCalciumDecayTau, 1 / (1 - (preCalciumRiseTau / preCalciumDecayTau))) -
                                   std::pow(preCalciumRiseTau / preCalciumDecayTau, 1 / ((preCalciumDecayTau / preCalciumRiseTau) - 1)))),
                                 -1);
  this->constants.postCalciumFluxFactor =
      postspikeCalcium * std::pow((((1 / postCalciumDecayTau) - (1 / postCalciumRiseTau)) *
                                   (std::pow(postCalciumRiseTau / postCalciumDecayTau, 1 / (1 - (postCalciumRiseTau / postCalciumDecayTau))) -
                                    std::pow(postCalciumRiseTau / postCalciumDecayTau, 1 / ((postCalciumDecayTau / postCalciumRiseTau) - 1)))),
                                  -1);
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
      throw "calciumBuffering was not consistent in plasticity model parameters.";
    } else if (parameterName.find("initialWeight") != std::string::npos && this->constants.initialWeight != std::stod(parameterValues.at(0))) {
      throw "initialWeight was not consistent in plasticity model parameters.";
    } else if (parameterName.find("availResourcesRatio") != std::string::npos && this->availResourcesRatio != std::stod(parameterValues.at(0))) {
      throw "availResourcesRatio was not consistent in plasticity model parameters.";
    } else if (parameterName.find("prespikeCalcium") != std::string::npos && this->prespikeCalcium != std::stod(parameterValues.at(0))) {
      throw "prespikeCalcium was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postspikeCalcium") != std::string::npos && this->postspikeCalcium != std::stod(parameterValues.at(0))) {
      throw "postspikeCalcium was not consistent in plasticity model parameters.";
    } else if (parameterName.find("preCalciumRiseTau") != std::string::npos &&
               this->preCalciumRiseTau != (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("preCalciumDecayTau") != std::string::npos &&
               this->preCalciumDecayTau != (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postCalciumRiseTau") != std::string::npos &&
               this->postCalciumRiseTau != (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("postCalciumDecayTau") != std::string::npos &&
               this->postCalciumDecayTau != (std::stod(parameterValues.at(0)) * infoGlobal->dtTimestep)) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("nonlinearNMDAFactor") != std::string::npos &&
               this->constants.nonlinearFactorNMDA != (std::stod(parameterValues.at(0)))) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
    } else if (parameterName.find("resourceConversionFct") != std::string::npos &&
               this->constants.resourceConversionFct != (std::stod(parameterValues.at(0)))) {
      throw "#mV*uM^-1 Conversion factor from nM to mV";

    } else if (parameterName.find("calciumBasal") != std::string::npos && this->calciumBasal != (std::stod(parameterValues.at(0)))) {
      throw "preCalciumDelay was not consistent in plasticity model parameters.";
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
  wParameterFile << neuronIdentificator << "preCalciumRiseTau\t\t" << (static_cast<double>(this->preCalciumRiseTau) / infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "preCalciumDecayTau\t\t" << (static_cast<double>(this->preCalciumDecayTau) / infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "postCalciumRiseTau\t\t" << (static_cast<double>(this->postCalciumRiseTau) / infoGlobal->dtTimestep)
                 << " #secs.";
  wParameterFile << "\t"
                 << "#Delay of the prespike calcium influx\n";
  wParameterFile << neuronIdentificator << "postCalciumDecayTau\t\t" << (static_cast<double>(this->postCalciumDecayTau) / infoGlobal->dtTimestep)
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
  wParameterFile << neuronIdentificator << "resourceConversionFct\t\t" << (constants.resourceConversionFct) << " #mV s^{-1} uM^{-1}";
  wParameterFile << "\t"
                 << "#Conversion from resource molar units (uM) to dmV/sec\n";
}

int MACRbPModel::CreateBranch(std::vector<int> anteriorBranches) {
  int branchId{this->GenerateBranchId()};
  if (!anteriorBranches.empty()) {
    this->caDiffBranches.emplace_back(MACRbPBranch(anteriorBranches, this->synapticGap, this->branchLength, branchId,
                                                   constants)); // This vector should be sorted by ID by default (tested).
    this->branches.push_back(static_cast<Branch *>(&this->caDiffBranches.back()));
  } else {
    int branchId{this->GenerateBranchId()};
    this->caDiffBranches.emplace_back(
        MACRbPBranch(this->synapticGap, this->branchLength, branchId, constants)); // This vector should be sorted by ID by default (tested).
    this->branches.push_back(static_cast<Branch *>(&this->caDiffBranches.back()));
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
  for (MACRbPBranch &branch : caDiffBranches) {
    branch.Advect();
    // branch.Advect(TStepModded);
  }
  // TStepModded++;//This is for the postspike calcium influx next timestep
}

void MACRbPModel::RecordPostSpike() {
  // TStepModded=infoGlobal->timeStep%preSpikeDelaySteps;
  this->totalPostSpikes++;
  // this->postSpiked = true;
  for (MACRbPBranch &branch : caDiffBranches) {
    branch.PostSpikeCalciumFlux();
  }
  // for (CaDiffusionBranch& branch: caDiffBranches){
  // 	std::transform(PAR_UNSEQ,branch.waitingMatrix.at(TStepModded).begin(), branch.waitingMatrix.at(TStepModded).end(),
  // branch.waitingMatrix.at(TStepModded).begin(),std::bind(std::plus<double>(),std::placeholders::_1, postspikeCalcium));
  // }
}

void MACRbPModel::RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) {
  // CaResSynapseSpine& synapseSpine = *caResSpines.at(spikedSpineId);
  // CaDiffusionBranch&  branch       = caDiffBranches.at(synapseSpine.branchId);
  // int  branchSpinePosition{synapseSpine.branchPositionId};
  static_cast<CaRsSpinePtr>(spinePtr)->preTransientIncrease++;
  // branch.waitingMatrix.at(TStepInput).at(branchSpinePosition)+=prespikeCalcium;
  this->totalPreSpikes++;
}

void MACRbPModel::PostConnectSetUp() {
  for (MACRbPBranch &branch : caDiffBranches) {
    branch.PostConnectSetUp(branchedSpineData);
  }
}

BaseSpinePtr MACRbPModel::AllocateNewSynapse(BranchTargeting &branchTarget) {
  int branch{AllocateBranch(branchTarget)};
  int position{PopSynapseSlotFromBranch(branchTarget)};
  // caDiffBranches.at(branch).CaDiffSpines.push_back(CaResSynapseSpine(kinasesTotal, calcineurinTotal, constants.initialWeight));
  CaRsSpinePtr newSpine = &caDiffBranches.at(branch).CaResSpines.at(position);

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
  this->caResSpines.push_back(newSpine);

  return this->baseSpineData.back();
}

std::vector<double> MACRbPModel::GetOverallSynapticProfile() const {
  return std::vector<double>();
}

std::string MACRbPModel::GetOverallSynapticProfileHeaderInfo() const {
  return std::string("{None}");
}
