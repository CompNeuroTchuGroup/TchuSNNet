//
// Created by Saif Ahmed on 15.07.21.
//
#include <string>
#include "MonoDendriteSTDPTazerartRelative.hpp"

MonoDendriteSTDPTazerartRelative::MonoDendriteSTDPTazerartRelative(GlobalSimInfo *info) : MonoDendriteSTDP(info) {
}

void MonoDendriteSTDPTazerartRelative::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    //Missing SP
    MonoDendriteSTDP::SaveParameters(stream, neuronPreId);
}

void MonoDendriteSTDPTazerartRelative::LoadParameters(std::vector<std::string> *input) {
    MonoDendriteSTDP::LoadParameters(input);

    std::string name;
    std::vector<std::string> values;

    bool muLTPInitialized = false,
            sigmaLTPInitialized = false,
            alphaInitialized = false,
            muLTDInitialized = false,
            sigmaLTDInitialized = false,
            betaInitialized = false;

    for (auto & it : *input) {
        SplitString(&it, &name, &values);

        if (name.find("mu_ltp") != std::string::npos) {
            this->muLTP = std::stod(values.at(0));
            muLTPInitialized = true;
        } else if (name.find("sigma_ltp") != std::string::npos) {
            this->sigmaLTP = std::stod(values.at(0));
            sigmaLTPInitialized = true;
        } else if (name.find("alpha_ltp") != std::string::npos) {
            this->alpha = std::stod(values.at(0));
            alphaInitialized = true;
        } else if (name.find("mu_ltd") != std::string::npos) {
            this->muLTD = std::stod(values.at(0));
            muLTDInitialized = true;
        } else if (name.find("sigma_ltd") != std::string::npos) {
            this->sigmaLTD = std::stod(values.at(0));
            sigmaLTDInitialized = true;
        } else if (name.find("beta_ltd") != std::string::npos) {
            this->beta = std::stod(values.at(0));
            betaInitialized = true;
        }
    }

    assertm(muLTPInitialized == true, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(sigmaLTPInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
    assertm(alphaInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");

    assertm(muLTDInitialized == true, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(sigmaLTDInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
    assertm(betaInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
}

void MonoDendriteSTDPTazerartRelative::updateLTP(unsigned long synId) {
    CoopSynapseSpine* syn = this->synapseDataCoop[synId].get();
//    this->weightsSum -= this->synapseDataCoop[synId]->weight;
    this->synapseDataCoop[synId]->AddToWeight(this->synapseDataCoop[synId]->GetWeight() * this->preFactorLTP * this->aLTP(syn->GetTheta()) * this->gLTP(this->lastPostSpikeTime - syn->GetLastSpike()));
//    this->synapseDataCoop[synId]->weight = std::min(2.0, this->synapseDataCoop[synId]->weight);
//    this->weightsSum += this->synapseDataCoop[synId]->weight;
}

void MonoDendriteSTDPTazerartRelative::updateLTD(unsigned long synId) {
    CoopSynapseSpine* syn = this->synapseDataCoop[synId].get();
//    this->weightsSum -= this->synapseDataCoop[synId]->weight;
    this->synapseDataCoop[synId]->AddToWeight(this->synapseDataCoop[synId]->GetWeight() * this->preFactorLTD * this->aLTD(syn->GetTheta()) * this->gLTD(syn->GetLastSpike() - this->lastPostSpikeTime));
//    this->synapseDataCoop[synId]->weight = std::max(0.0, this->synapseDataCoop[synId]->weight);
//    this->weightsSum += this->synapseDataCoop[synId]->weight;
}

double MonoDendriteSTDPTazerartRelative::gLTP(double deltaT) const {
    double denominator = 2 * this->sigmaLTP * this->sigmaLTP;
    double numerator = this->muLTP - deltaT;
    numerator *= numerator;
    return exp(-numerator / denominator);
}

double MonoDendriteSTDPTazerartRelative::gLTD(double deltaT) const {
    double denominator = 2 * this->sigmaLTD * this->sigmaLTD;
    double numerator = this->muLTD - deltaT;
    numerator *= numerator;
    return exp(-numerator / denominator);
}

double MonoDendriteSTDPTazerartRelative::aLTP(double theta) const {
    return 2.0 - exp(-this->alpha * theta);
}

double MonoDendriteSTDPTazerartRelative::aLTD(double theta) const {
    return -exp(-this->beta * theta);
}

const std::string MonoDendriteSTDPTazerartRelative::GetType() {
    return str_MonoDendriteSTDPTazerartRelative;
}

double MonoDendriteSTDPTazerartRelative::getTimingEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const {
    if (synA == synB) {
        return 0.0;
    }
    if (synA->GetLastSpike() < 0 || synB->GetLastSpike() < 0) {
        return 0.0;
    }
    return exp(-abs(synA->GetLastSpike() - synB->GetLastSpike()) / this->tauDelay);
}

double MonoDendriteSTDPTazerartRelative::getDistanceEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const {
    if (synA == synB) {
        return 0;
    }
    return exp(-abs(synA->GetDistToSoma() - synB->GetDistToSoma()) / this->lambdaDist);
}
