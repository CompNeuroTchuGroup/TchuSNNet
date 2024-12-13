//
// Created by Saif Ahmed on 15.07.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#include "MonoDendriteSTDPTazerart.hpp"

MonoDendriteSTDPTazerart::MonoDendriteSTDPTazerart(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters)
    : MonoDendriteSTDP(infoGlobal, morphologyParameters) {
}

void MonoDendriteSTDPTazerart::SaveParameters(std::ofstream &wParameterFile, std::string neuronIdentificator) const {
  MonoDendriteSTDP::SaveParameters(wParameterFile, neuronIdentificator);
  wParameterFile << "#Tazerart's exclusive here.\n";
  wParameterFile << neuronIdentificator << "muLtp\t\t\t\t" << std::to_string(this->muLTP);
  wParameterFile << "\t"
                 << "#Time interval at which the LTP effect is maximized.\n";

  wParameterFile << neuronIdentificator << "sigmaLtp\t\t\t\t" << std::to_string(this->sigmaLTP);
  wParameterFile << "\t"
                 << "#used in gLTP, spread of the effect distribution.\n";

  wParameterFile << neuronIdentificator << "alphaLtp\t\t\t\t" << std::to_string(this->alpha);
  wParameterFile << "\t"
                 << "#used in aLTP.\n";

  wParameterFile << neuronIdentificator << "muLtd\t\t\t\t" << std::to_string(this->muLTD);
  wParameterFile << "\t"
                 << "#(negative) interval of time at which the ltd effect is maximized.\n";

  wParameterFile << neuronIdentificator << "sigmaLtd\t\t\t\t" << std::to_string(this->sigmaLTD);
  wParameterFile << "\t"
                 << "#used in gLTD, spread of the effect distribution.\n";

  wParameterFile << neuronIdentificator << "betaLtd\t\t\t\t" << std::to_string(this->beta);
  wParameterFile << "\t"
                 << "#used in aLTD.\n";
}

void MonoDendriteSTDPTazerart::LoadParameters(const std::vector<FileEntry> &morphologyParameters) {
  MonoDendriteSTDP::LoadParameters(morphologyParameters);

  bool muLTPInitialized = false, sigmaLTPInitialized = false, muLTDInitialized = false, sigmaLTDInitialized = false;

  for (auto &[parameterName, parameterValues] : morphologyParameters) {
    if (parameterName.find("muLtp") != std::string::npos) {
      this->muLTP      = std::stod(parameterValues.at(0));
      muLTPInitialized = true;
    } else if (parameterName.find("sigmaLtp") != std::string::npos) {
      this->sigmaLTP      = std::stod(parameterValues.at(0));
      sigmaLTPInitialized = true;
    } else if (parameterName.find("muLtd") != std::string::npos) {
      this->muLTD      = std::stod(parameterValues.at(0));
      muLTDInitialized = true;
    } else if (parameterName.find("sigmaLtd") != std::string::npos) {
      this->sigmaLTD      = std::stod(parameterValues.at(0));
      sigmaLTDInitialized = true;
    }
  }

  if (!muLTPInitialized) {
    throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
  } else if (!sigmaLTPInitialized) {
    throw "Using heterosynaptic synapses without specifying synapticGap is not allowed.";
  } else if (!muLTDInitialized) {
    throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
  } else if (!sigmaLTDInitialized) {
    throw "Using heterosynaptic synapses without specifying synapticGap is not allowed.";
  }
}

void MonoDendriteSTDPTazerart::CheckParameters(const std::vector<FileEntry> &parameters) {
  for (auto &[parameterName, parameterValues] : parameters) {
    if (parameterName.find("muLtp") != std::string::npos) {
      if (!(this->muLTP == std::stod(parameterValues.at(0)))) {
        throw "dendriticLength was not consistent in plasticity parameters";
      }
    } else if (parameterName.find("sigmaLtp") != std::string::npos) {
      if (!(this->sigmaLTP == std::stod(parameterValues.at(0)))) {
        throw "sigmaLtp was not consistent in plasticity parameters";
      }
    } else if (parameterName.find("muLtd") != std::string::npos) {
      if (!(this->muLTD == std::stod(parameterValues.at(0)))) {
        throw "muLtd was not consistent in plasticity parameters";
      }
    } else if (parameterName.find("sigmaLtd") != std::string::npos) {
      if (!(this->sigmaLTD == std::stod(parameterValues.at(0)))) {
        throw "sigmaLtd was not consistent in plasticity parameters";
      }
    }
  }
}

void MonoDendriteSTDPTazerart::UpdateLTP(signed long spineID) {
  CoopSynapseSpine *spine = this->spineDataCoop.at(spineID).get();
  //    this->weightsSum -= this->synapseDataCoop.at(synId)->weight;
  double change = this->preFactorLTP * this->aLTP(spine->GetTheta()) * this->gLTP(this->lastPostSpikeTime - spine->GetLastSpike());
  //    std::cout << synId << " : "  << change << std::endl;
  this->spineDataCoop.at(spineID)->weight += (change);

  //    if (synId == 0) {
  //        std::cout << synId << " -- "<< this->lastPostSpikeTime << ", " << syn->lastSpike << " : "  << change << std::endl;
  //    }
  //    this->synapseDataCoop.at(synId)->weight = std::min(2.0, this->synapseDataCoop.at(synId)->weight);
  //    this->weightsSum += this->synapseDataCoop.at(synId)->weight;
  //    this->ltpList.push_back(change);
  //    for (signed long i = 0; i < this->ltpList.size(); i++) {
  //        std::cout << this->ltpList.at(i) << ", ";
  //    }
  //    std::cout << std::endl;
}

void MonoDendriteSTDPTazerart::UpdateLTD(signed long spineID) {
  CoopSynapseSpine *spine = this->spineDataCoop.at(spineID).get();
  //    this->weightsSum -= this->synapseDataCoop.at(synId)->weight;
  double change = -this->preFactorLTD * this->aLTD(spine->GetTheta()) * this->gLTD(spine->GetLastSpike() - this->lastPostSpikeTime);
  //    std::cout << synId << " : "  << change << std::endl;
  this->spineDataCoop.at(spineID)->weight += (change);

  //    if (synId == 0) {
  //        std::cout << synId << " -- " << syn->lastSpike << " , " << this->lastPostSpikeTime << " : "  << change << std::endl;
  //    }
  //    this->synapseDataCoop.at(synId)->weight = std::max(0.0, this->synapseDataCoop.at(synId)->weight);
  //    this->weightsSum += this->synapseDataCoop.at(synId)->weight;
  //    this->ltdList.push_back(change);
  //    for (signed long i = 0; i < this->ltdList.size(); i++) {
  //        std::cout << this->ltdList.at(i) << ", ";
  //    }
  //    std::cout << std::endl;
}

double MonoDendriteSTDPTazerart::gLTP(double deltaT) const {
  if (deltaT < 0.0) {
    return 0.0;
  } else {
    return exp(-std::pow(this->muLTP - deltaT, 2) / (2 * std::pow(this->sigmaLTP, 2)));
  }
}

double MonoDendriteSTDPTazerart::gLTD(double deltaT) const {
  if (deltaT <= 0.0) {
    return 0.0;
  } else {
    return exp(-std::pow(this->muLTD - deltaT, 2) / (2 * std::pow(this->sigmaLTD, 2)));
  }
}