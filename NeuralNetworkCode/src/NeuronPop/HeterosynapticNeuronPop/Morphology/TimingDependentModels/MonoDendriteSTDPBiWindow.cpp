//
// Created by Saif Ahmed on 15.07.21.
//
#include <string>
#include "MonoDendriteSTDPBiWindow.hpp"

MonoDendriteSTDPBiWindow::MonoDendriteSTDPBiWindow(GlobalSimInfo *info) : MonoDendriteSTDP(info) {
}

void MonoDendriteSTDPBiWindow::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    MonoDendriteSTDP::SaveParameters(stream, neuronPreId);

    *stream << neuronPreId<<"_tau_ltp\t\t\t"<<std::to_string(this->tauLTP);
    *stream << "\t"<<"#Decay constant of the temporal effect in LTP.\n";
    *stream << neuronPreId<<"_tau_ltd\t\t\t"<<std::to_string(this->tauLTD);
    *stream << "\t"<<"#Decay constant of the temporal effect in LTD.\n";
    *stream << neuronPreId<<"_alpha_ltp\t\t\t"<<std::to_string(this->alpha);
    *stream << "\t"<<"#Cooperativity decay for LTP.\n";
    *stream << neuronPreId<<"_beta_ltd\t\t\t"<<std::to_string(this->beta);
    *stream << "\t"<<"#Cooperativity decay for LTD.\n";

}

void MonoDendriteSTDPBiWindow::LoadParameters(std::vector<std::string> *input) {
    MonoDendriteSTDP::LoadParameters(input);

    std::string name;
    std::vector<std::string> values;

    bool tauLTPInitialized = false, tauLTDInitialized = false, alphaInitialized = false, betaInitialized = false;

    for (auto & it : *input) {
        SplitString(&it, &name, &values);

        if (name.find("tau_ltp") != std::string::npos) {
            this->tauLTP = std::stod(values.at(0));
            tauLTPInitialized = true;
        } else if (name.find("tau_ltd") != std::string::npos) {
            this->tauLTD = std::stod(values.at(0));
            tauLTDInitialized = true;
        }  else if (name.find("alpha_ltp") != std::string::npos) {
            this->alpha = std::stod(values.at(0));
            alphaInitialized = true;
        } else if (name.find("beta_ltd") != std::string::npos) {
            this->beta = std::stod(values.at(0));
            betaInitialized = true;
        }
    }

    assertm(tauLTPInitialized == true, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(alphaInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");

    assertm(tauLTDInitialized == true, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(betaInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
}

void MonoDendriteSTDPBiWindow::updateLTP(unsigned long synId) {
    CoopSynapseSpine* syn = this->synapseDataCoop[synId].get();
//    this->weightsSum -= this->synapseDataCoop[synId]->weight;
//    double time = this->lastPostSpikeTime - syn->GetLastSpike();
    double change = this->preFactorLTP * this->aLTP(syn->GetTheta()) * this->gLTP(this->lastPostSpikeTime - syn->GetLastSpike());
    this->synapseDataCoop[synId]->AddToWeight(change);

//    this->weight_changes.emplace_back(synId, change);

//    if (synId == 0) {
//        std::cout << change << std::endl;
//    }

//    this->synapseDataCoop[synId]->weight = std::min(2.0, this->synapseDataCoop[synId]->weight);
//    this->weightsSum += this->synapseDataCoop[synId]->weight;

}

void MonoDendriteSTDPBiWindow::updateLTD(unsigned long synId) {
    CoopSynapseSpine* syn = this->synapseDataCoop[synId].get();
//    this->weightsSum -= this->synapseDataCoop[synId]->weight;

//    double time = this->lastPostSpikeTime - syn->lastSpike;
    double change  = this->preFactorLTD * this->aLTD(syn->GetTheta()) * this->gLTD(this->lastPostSpikeTime - syn->GetLastSpike());
    this->synapseDataCoop[synId]->AddToWeight(change);

//    this->weight_changes.emplace_back(synId, change);


//    std::cout << change << std::endl;

//    this->synapseDataCoop[synId]->weight = std::max(0.0, this->synapseDataCoop[synId]->weight);
//    this->weightsSum += this->synapseDataCoop[synId]->weight;
}

double MonoDendriteSTDPBiWindow::gLTP(double deltaT) const {
    if (deltaT <= 0.0) return 0.0;
    return exp(-deltaT / this->tauLTP);
}

double MonoDendriteSTDPBiWindow::gLTD(double deltaT) const {
    if (deltaT >= 0.0) return 0.0;
    return -exp(deltaT / this->tauLTD);
}

double MonoDendriteSTDPBiWindow::aLTP(double theta) const {//coop
    return base_ltp +  incr_ltp * (1 - exp(-this->alpha * theta));
}

double MonoDendriteSTDPBiWindow::aLTD(double theta) const {//coop
    return (base_ltd - decr_ltd * (1 - exp(-this->beta * theta)));
}

const std::string MonoDendriteSTDPBiWindow::GetType() {
    return str_MonoDendriteSTDPBiWindow;
}

double MonoDendriteSTDPBiWindow::getTimingEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const {
    if (synA == synB) {
        return 0.0;
    }
    if (synA->GetLastSpike() < 0 || synB->GetLastSpike() < 0) {
        return 0.0;
    }
    return exp(-std::abs(synA->GetLastSpike() - synB->GetLastSpike()) / this->tauDelay);
}

double MonoDendriteSTDPBiWindow::getDistanceEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const {
    if (synA == synB) {
        return 0;
    }
    return exp(-std::abs(synA->GetDistToSoma() - synB->GetDistToSoma()) / this->lambdaDist);
}
