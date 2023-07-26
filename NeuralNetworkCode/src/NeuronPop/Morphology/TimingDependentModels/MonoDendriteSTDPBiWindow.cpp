//
// Created by Saif Ahmed on 15.07.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#include "MonoDendriteSTDPBiWindow.hpp"
#include <string>


MonoDendriteSTDPBiWindow::MonoDendriteSTDPBiWindow(GlobalSimInfo* infoGlobal) : MonoDendriteSTDP(infoGlobal) {
}

void MonoDendriteSTDPBiWindow::SaveParameters(std::ofstream& wParameterFile, std::string neuronIDString) const {
    MonoDendriteSTDP::SaveParameters(wParameterFile, neuronIDString);

    wParameterFile << neuronIDString<<"tauLtp\t\t\t\t"<<std::to_string(this->tauLTP);
    wParameterFile << "\t"<<"#Decay constant of the temporal effect in LTP.\n";
    wParameterFile << neuronIDString<<"tauLtd\t\t\t\t"<<std::to_string(this->tauLTD);
    wParameterFile << "\t"<<"#Decay constant of the temporal effect in LTD.\n";
}

void MonoDendriteSTDPBiWindow::LoadParameters(const std::vector<FileEntry>& morphologyParameters) {//const
    MonoDendriteSTDP::LoadParameters(morphologyParameters);

    bool tauLTPInitialized = false, tauLTDInitialized = false;

    for (auto& [parameterName, parameterValues] : morphologyParameters) {
        if (parameterName.find("tauLtp") != std::string::npos) {
            this->tauLTP = std::stod(parameterValues.at(0));
            tauLTPInitialized = true;
        } else if (parameterName.find("tauLtd") != std::string::npos) {
            this->tauLTD = std::stod(parameterValues.at(0));
            tauLTDInitialized = true;
        }
    }
    if (!tauLTDInitialized){
        throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
    } else if (!tauLTPInitialized){
        throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
    }
}

void MonoDendriteSTDPBiWindow::CheckParameters(const std::vector<FileEntry> &parameters) {
    MonoDendriteSTDP::CheckParameters(parameters);
    for (auto& [parameterName, parameterValues] : parameters) {
        if (parameterName.find("tauLtp") != std::string::npos) {
            if(!(this->tauLTP == std::stod(parameterValues.at(0)))){
                throw "tauLtp was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("tauLtd") != std::string::npos) {
            if(!(this->tauLTD == std::stod(parameterValues.at(0)))){
                throw "tauLtd was not consistent in plasticity parameters";
            }
        }
    }
}

void MonoDendriteSTDPBiWindow::UpdateLTP(signed long spineID) {
    CoopSynapseSpine* spine = this->spineDataCoop.at(spineID).get();
//    this->weightsSum -= this->synapseDataCoop.at(synId)->weight;
//    double time = this->lastPostSpikeTime - syn->GetLastSpike();
    double change = this->preFactorLTP * this->aLTP(spine->GetTheta()) * this->gLTP(this->lastPostSpikeTime - spine->GetLastSpike());
    this->spineDataCoop.at(spineID)->AddToWeight(change);

//    this->weight_changes.emplace_back(synId, change);

//    if (synId == 0) {
//        std::cout << change << std::endl;
//    }

//    this->synapseDataCoop.at(synId)->weight = std::min(2.0, this->synapseDataCoop.at(synId)->weight);
//    this->weightsSum += this->synapseDataCoop.at(synId)->weight;

}

void MonoDendriteSTDPBiWindow::UpdateLTD(signed long spineID) {
    CoopSynapseSpine* spine = this->spineDataCoop.at(spineID).get();
//    this->weightsSum -= this->synapseDataCoop.at(synId)->weight;

//    double time = this->lastPostSpikeTime - syn->lastSpike;
    double change  = -this->preFactorLTD * this->aLTD(spine->GetTheta()) * this->gLTD(this->lastPostSpikeTime - spine->GetLastSpike());
    this->spineDataCoop.at(spineID)->AddToWeight(change);

//    this->weight_changes.emplace_back(synId, change);


//    std::cout << change << std::endl;

//    this->synapseDataCoop.at(synId)->weight = std::max(0.0, this->synapseDataCoop.at(synId)->weight);
//    this->weightsSum += this->synapseDataCoop.at(synId)->weight;
}

double MonoDendriteSTDPBiWindow::gLTP(double deltaT) const {
    if (deltaT <= 0.0) {return 0.0;}
    return exp(-deltaT / this->tauLTP);
}

double MonoDendriteSTDPBiWindow::gLTD(double deltaT) const {
    if (deltaT >= 0.0) {return 0.0;}
    return exp(deltaT / this->tauLTD);
}
