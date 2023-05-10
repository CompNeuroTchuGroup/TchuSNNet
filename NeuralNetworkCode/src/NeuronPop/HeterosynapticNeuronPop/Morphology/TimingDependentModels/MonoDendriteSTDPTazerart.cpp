//
// Created by Saif Ahmed on 15.07.21.
//
#include <string>
#include "MonoDendriteSTDPTazerart.hpp"

MonoDendriteSTDPTazerart::MonoDendriteSTDPTazerart(GlobalSimInfo *info) : MonoDendriteSTDP(info) {
}

void MonoDendriteSTDPTazerart::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    MonoDendriteSTDP::SaveParameters(stream, neuronPreId);
    *stream <<"#Tazerart's exclusive here.\n";
    *stream << neuronPreId<<"_mu_ltp\t\t\t\t\t"<<std::to_string(this->muLTP);
    *stream << "\t"<<"#Time interval at which the LTP effect is maximized.\n";

    *stream << neuronPreId<<"_sigma_ltp\t\t\t\t"<<std::to_string(this->sigmaLTP);
    *stream << "\t"<<"#used in gLTP, spread of the effect distribution.\n";

    *stream << neuronPreId<<"_alpha_ltp\t\t\t\t"<<std::to_string(this->alpha);
    *stream << "\t"<<"#used in aLTP.\n";

    *stream << neuronPreId<<"_mu_ltd\t\t\t\t\t"<<std::to_string(this->muLTD);
    *stream << "\t"<<"#(negative) interval of time at which the ltd effect is maximized.\n";

    *stream << neuronPreId<<"_sigma_ltd\t\t\t\t"<<std::to_string(this->sigmaLTD);
    *stream << "\t"<<"#used in gLTD, spread of the effect distribution.\n";

    *stream << neuronPreId<<"_beta_ltd\t\t\t\t"<<std::to_string(this->beta);
    *stream << "\t"<<"#used in aLTD.\n";
}

void MonoDendriteSTDPTazerart::LoadParameters(std::vector<std::string> *input) {
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

void MonoDendriteSTDPTazerart::updateLTP(unsigned long synId) {
    CoopSynapseSpine* syn = this->synapseDataCoop[synId].get();
//    this->weightsSum -= this->synapseDataCoop[synId]->weight;
    double change = this->preFactorLTP * this->aLTP(syn->GetTheta()) * this->gLTP(this->lastPostSpikeTime - syn->GetLastSpike());
//    std::cout << synId << " : "  << change << std::endl;
    this->synapseDataCoop[synId]->AddToWeight(change);

//    if (synId == 0) {
//        std::cout << synId << " -- "<< this->lastPostSpikeTime << ", " << syn->lastSpike << " : "  << change << std::endl;
//    }
//    this->synapseDataCoop[synId]->weight = std::min(2.0, this->synapseDataCoop[synId]->weight);
//    this->weightsSum += this->synapseDataCoop[synId]->weight;
//    this->ltpList.push_back(change);
//    for (unsigned long i = 0; i < this->ltpList.size(); i++) {
//        std::cout << this->ltpList.at(i) << ", ";
//    }
//    std::cout << std::endl;
}


void MonoDendriteSTDPTazerart::updateLTD(unsigned long synId) {
    CoopSynapseSpine* syn = this->synapseDataCoop[synId].get();
//    this->weightsSum -= this->synapseDataCoop[synId]->weight;
    double change = this->preFactorLTD * this->aLTD(syn->GetTheta()) * this->gLTD(syn->GetLastSpike() - this->lastPostSpikeTime);
//    std::cout << synId << " : "  << change << std::endl;
    this->synapseDataCoop[synId]->AddToWeight(change);

//    if (synId == 0) {
//        std::cout << synId << " -- " << syn->lastSpike << " , " << this->lastPostSpikeTime << " : "  << change << std::endl;
//    }
//    this->synapseDataCoop[synId]->weight = std::max(0.0, this->synapseDataCoop[synId]->weight);
//    this->weightsSum += this->synapseDataCoop[synId]->weight;
//    this->ltdList.push_back(change);
//    for (unsigned long i = 0; i < this->ltdList.size(); i++) {
//        std::cout << this->ltdList.at(i) << ", ";
//    }
//    std::cout << std::endl;
}

double MonoDendriteSTDPTazerart::gLTP(double deltaT) const {
    if (deltaT < 0.0) {
        return 0.0;
    }
    double denominator = 2 * this->sigmaLTP * this->sigmaLTP;
    double numerator = this->muLTP - deltaT;
    numerator *= numerator;
    return exp(-numerator / denominator);
}

double MonoDendriteSTDPTazerart::gLTD(double deltaT) const {
    if (deltaT <= 0.0) {
        return 0.0;
    }
    double denominator = 2 * this->sigmaLTD * this->sigmaLTD;
    double numerator = this->muLTD - deltaT;
    numerator *= numerator;
    return exp(-numerator / denominator);
}

double MonoDendriteSTDPTazerart::aLTP(double theta) const {
    return base_ltp + incr_ltp * (1 - exp(-this->alpha * theta));
}

double MonoDendriteSTDPTazerart::aLTD(double theta) const {
    return - (base_ltd - decr_ltd * (1 - exp(-this->beta * theta)));
}

const std::string MonoDendriteSTDPTazerart::GetType() {
    return str_MonoDendriteSTDPTazerart;
}

double MonoDendriteSTDPTazerart::getTimingEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const {
    if (synA == synB) {
        return 0.0;
    }
    if (synA->GetLastSpike() < 0 || synB->GetLastSpike() < 0) {
        return 0.0;
    }
    return exp(-std::abs(synA->GetLastSpike() - synB->GetLastSpike()) / this->tauDelay);
}

double MonoDendriteSTDPTazerart::getDistanceEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const {
    if (synA == synB) {
        return 0;
    }
    return exp(-std::abs(synA->GetDistToSoma() - synB->GetDistToSoma()) / this->lambdaDist);
}
