//
// Created by Saif Ahmed on 05.08.21.
//

#include "HeteroPoissonNeuronPop.hpp"

HeteroPoissonNeuronPop::HeteroPoissonNeuronPop(GlobalSimInfo *info, int id): HeteroNeuronPop(info, id)    {
    r_target = 0; seed = 2;
    generator = std::mt19937(seed);
    uniformDistribution = std::uniform_real_distribution<double>(0.0,1.0);
}

void HeteroPoissonNeuronPop::advect(std::vector<double> * synaptic_dV)
{
    // double dt           = info->dt;

    ClearSpiker();

    //#pragma omp parallel for
    for(unsigned long i = 0 ; i < noNeurons; i++)
    {
        // set target rate
        if (inputDependant){
            r_target = synaptic_dV->at(i);
            lambda = r_target;
        }

        //check spiking
        if (uniformDistribution(generator) < lambda)
            spiker.push_back(i);
    }

    for (auto neuron: this->spiker) {
        this->morphology[neuron]->RecordPostSpike();
    }

    for (unsigned long morphId =  0; morphId < morphology.size(); ++morphId) {
        this->morphology[morphId]->advect();
    }


}

void HeteroPoissonNeuronPop::LoadParameters(std::vector<std::string> *input){

    HeteroNeuronPop::LoadParameters(input);

    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);
        if(name.find("r_target") != std::string::npos){
            r_target = std::stoi(values.at(0));
            inputDependant = false;
            lambda 	= r_target*this->info->dt;
        } else if (name.find("seedPoisson") != std::string::npos){
            seed = std::stoi(values.at(0));
        }
    }

    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        seed = distribution(info->globalGenerator);
        generator = std::mt19937(seed);
    }

}


void HeteroPoissonNeuronPop::SaveParameters(std::ofstream * stream){

    std::string id = "neurons_" + std::to_string(GetId());

    //std::cout<<"printing neuron id: "<<id<<"\n";

    //NeuronPop::SaveParameters(stream);
    *stream <<  "#***********************************************\n";
    *stream <<  id + "_noNeurons                   " << noNeurons << "\n";
    *stream <<  id + "_type                        " << GetType() << "\n";
    if (info->globalSeed == -1) {
        *stream << id + "_seedPoisson                 " << std::to_string(seed) << "\n";
    }
    if (!inputDependant){
        *stream <<  id + "_r_target                   " << std::to_string(r_target)  << "\n";
    }
    *stream <<  "#\t\tHeteroPoisson neuron: produces Poisson spiking with rate r_target (defined under stimulus) \n";

}

std::string HeteroPoissonNeuronPop::GetType() {
    return str_HeteroPoissonNeuron;
}
