//
// Created by Saif Ahmed on 15.07.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#include <string>
#include "MonoDendriteSTDPTazerartRelative.hpp"

MonoDendriteSTDPTazerartRelative::MonoDendriteSTDPTazerartRelative(GlobalSimInfo* infoGlobal) : MonoDendriteSTDP(infoGlobal) {
}

void MonoDendriteSTDPTazerartRelative::SaveParameters(std::ofstream& wParameterFile, std::string neuronIdentificator) const {
    //Missing SP
    MonoDendriteSTDP::SaveParameters(wParameterFile, neuronIdentificator);
}

void MonoDendriteSTDPTazerartRelative::LoadParameters(const std::vector<FileEntry>& morphologyParameters) {
    MonoDendriteSTDP::LoadParameters(morphologyParameters);

    bool muLTPInitialized = false,
            sigmaLTPInitialized = false,
            muLTDInitialized = false,
            sigmaLTDInitialized = false;

    for (auto& [parameterName, parameterValues] : morphologyParameters) {
        if (parameterName.find("muLtp") != std::string::npos) {
            this->muLTP = std::stod(parameterValues.at(0));
            muLTPInitialized = true;
        } else if (parameterName.find("sigmaLtp") != std::string::npos) {
            this->sigmaLTP = std::stod(parameterValues.at(0));
            sigmaLTPInitialized = true;
        } else if (parameterName.find("muLtd") != std::string::npos) {
            this->muLTD = std::stod(parameterValues.at(0));
            muLTDInitialized = true;
        } else if (parameterName.find("sigmaLtd") != std::string::npos) {
            this->sigmaLTD = std::stod(parameterValues.at(0));
            sigmaLTDInitialized = true;
        }
    }
    if (!muLTPInitialized){
        throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
    } else if (!sigmaLTPInitialized){
        throw "Using heterosynaptic synapses without specifying synapticGap is not allowed.";
    } else if (!muLTDInitialized){
        throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
    } else if (!sigmaLTDInitialized){
        throw "Using heterosynaptic synapses without specifying synapticGap is not allowed.";
    }
}

void MonoDendriteSTDPTazerartRelative::CheckParameters(const std::vector<FileEntry> &parameters) {
    for (auto& [parameterName, parameterValues] : parameters) {
        if (parameterName.find("muLtp") != std::string::npos) {
            if(!(this->muLTP == std::stod(parameterValues.at(0)))){
                throw "dendriticLength was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("sigmaLtp") != std::string::npos) {
            if(!(this->sigmaLTP == std::stod(parameterValues.at(0)))){
                throw "sigmaLtp was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("muLtd") != std::string::npos) {
            if(!(this->muLTD == std::stod(parameterValues.at(0)))){
                throw "muLtd was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("sigmaLtd") != std::string::npos) {
            if(!(this->sigmaLTD == std::stod(parameterValues.at(0)))){
                throw "sigmaLtd was not consistent in plasticity parameters";
            }
        }
    }
}

void MonoDendriteSTDPTazerartRelative::UpdateLTP(signed long spineID) {
    CoopSynapseSpine* spine = this->spineDataCoop.at(spineID).get();
//    this->weightsSum -= this->synapseDataCoop.at(synId)->weight;
    this->spineDataCoop.at(spineID)->AddToWeight(this->spineDataCoop.at(spineID)->GetWeightUncoupled() * this->preFactorLTP * this->aLTP(spine->GetTheta()) * this->gLTP(this->lastPostSpikeTime - spine->GetLastSpike()));
//    this->synapseDataCoop.at(synId)->weight = std::min(2.0, this->synapseDataCoop.at(synId)->weight);
//    this->weightsSum += this->synapseDataCoop.at(synId)->weight;
}

void MonoDendriteSTDPTazerartRelative::UpdateLTD(signed long spineID) {
    CoopSynapseSpine* spine = this->spineDataCoop.at(spineID).get();
//    this->weightsSum -= this->synapseDataCoop.at(synId)->weight;
    this->spineDataCoop.at(spineID)->AddToWeight(-this->spineDataCoop.at(spineID)->GetWeightUncoupled() * this->preFactorLTD * this->aLTD(spine->GetTheta()) * this->gLTD(spine->GetLastSpike() - this->lastPostSpikeTime));
//    this->synapseDataCoop.at(synId)->weight = std::max(0.0, this->synapseDataCoop.at(synId)->weight);
//    this->weightsSum += this->synapseDataCoop.at(synId)->weight;
}

double MonoDendriteSTDPTazerartRelative::gLTP(double deltaT) const {
    return exp(-std::pow(this->muLTP - deltaT, 2) / (2 * std::pow(this->sigmaLTP, 2)));
}

double MonoDendriteSTDPTazerartRelative::gLTD(double deltaT) const {
    return exp(-std::pow(this->muLTD - deltaT, 2) / (2 * std::pow(this->sigmaLTD, 2)));
}

double MonoDendriteSTDPTazerartRelative::aLTP(double theta) const {
    return 2.0 - exp(-this->alpha * theta);
}

double MonoDendriteSTDPTazerartRelative::aLTD(double theta) const {
    return exp(-this->beta * theta);
}