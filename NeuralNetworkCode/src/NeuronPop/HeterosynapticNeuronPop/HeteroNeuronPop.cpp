//
// Created by Saif Ahmed on 03.08.21.
//

#include "HeteroNeuronPop.hpp"

HeteroNeuronPop::HeteroNeuronPop(GlobalSimInfo *info, int id): NeuronPop(info, id) {
}

void HeteroNeuronPop::SaveParameters(std::ofstream * stream) {
    NeuronPop::SaveParameters(stream);
    if (this->morphology.size()!=0){
        morphology.at(0)->SaveParameters(stream, "neurons_" + std::to_string(GetId()));
    }
}

void HeteroNeuronPop::LoadParameters(std::vector<std::string> *input) {

    NeuronPop::LoadParameters(input);

    std::string name;
    std::vector<std::string> values;


    for (auto & it : *input) {
        SplitString(&it,&name,&values);

        if(name.find("morphology_type") != std::string::npos) {
            if (values.at(0) == str_MonoDendriteSTDPTazerart) {
                for (unsigned long i = 0; i < this->noNeurons; i++) {
                    this->morphology.push_back(std::make_unique<MonoDendriteSTDPTazerart>(this->info));
                    this->morphology.back()->LoadParameters(input);
                }
            } else if (values.at(0) == str_MonoDendriteSTDPBiWindow) {
                for (unsigned long i = 0; i < this->noNeurons; i++) {
                    this->morphology.push_back(std::make_unique<MonoDendriteSTDPBiWindow>(this->info));
                    this->morphology.back()->LoadParameters(input);
                }
            } else if (values.at(0) == str_MonoDendriteSTDPTazerartRelative) {
                for (unsigned long i = 0; i < this->noNeurons; i++) {
                    this->morphology.push_back(std::make_unique<MonoDendriteSTDPTazerartRelative>(this->info));
                    this->morphology.back()->LoadParameters(input);
                }
            } else if (values.at(0) == str_BranchedResourceHeteroSTDP) {
                for (unsigned long i = 0; i < this->noNeurons; i++) {
                    this->morphology.push_back(std::make_unique<BranchedResourceHeteroSTDP>(this->info)); //Remove, will not  be used
                    this->morphology.back()->LoadParameters(input);
                }
            } else {
                throw;
            }
        }
    }

    if (this->morphology.at(0)->IsBranchedBool()){
        this->SetBranchedTrue();
    }
}

void HeteroNeuronPop::PostConnectSetUp()
{
    for (std::unique_ptr<Morphology>& singleMorphology : morphology){
        singleMorphology->PostConnectSetUp();
    }
}

BaseSpinePtr HeteroNeuronPop::AllocateNewSynapse(unsigned long neuronId, HeteroCurrentSynapse &synapse)
{
    return this->morphology[neuronId]->AllocateNewSynapse(synapse);
}

void HeteroNeuronPop::RecordExcitatorySynapticSpike(unsigned long neuronId, unsigned long synapseId)
{
        //This function is NOT DELAY COMPATIBLE
    this->morphology[neuronId]->RecordExcitatoryPreSpike(synapseId);
}

// std::vector<unsigned long> getSpikedSynapses(const HeteroNeuronPop& pop, unsigned long neuronId) {
//     return getSpikedSynapsesFromMorphology(*pop.morphology[neuronId].get());
// }

unsigned long HeteroNeuronPop::GetSynapseCount(unsigned long neuronId) {
    return this->morphology.at(neuronId)->GetSynapseCount();
}

double HeteroNeuronPop::GetWeight(unsigned long neuronId, unsigned long synapseId) {
    return this->morphology.at(neuronId)->GetWeight(synapseId);
}

std::valarray<double> HeteroNeuronPop::GetIndividualSynapticProfile(unsigned long neuronId, unsigned long synapseId) {
    return this->morphology.at(neuronId)->GetIndividualSynapticProfile(synapseId);
}

std::valarray<double> HeteroNeuronPop::GetOverallSynapticProfile(unsigned long neuronId) {
    return this->morphology.at(neuronId)->GetOverallSynapticProfile();
}
std::string HeteroNeuronPop::GetIndividualSynapticProfileHeaderInfo() const
{
    return morphology.at(0)->GetIndividualSynapticProfileHeaderInfo();
}
std::string HeteroNeuronPop::GetOverallSynapticProfileHeaderInfo() const
{
    return morphology.at(0)->GetOverallSynapticProfileHeaderInfo();
}
// STDP Analysis
// void HeteroNeuronPop::triggerStatOut(std::string dirPath) {
// //    this->morphology.at(0)->triggerStatOut(dirPath);
// }

// void HeteroNeuronPop::printThetasAndWeights() {
//     for (auto& m: this->morphology) {
//         m->printThetasAndWeights();
//     }
// }
