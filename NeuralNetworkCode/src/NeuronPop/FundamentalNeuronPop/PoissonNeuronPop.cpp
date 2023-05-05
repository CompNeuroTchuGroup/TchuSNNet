#include "PoissonNeuronPop.hpp"

PoissonNeuronPop::PoissonNeuronPop(GlobalSimInfo * info,int id):NeuronPop(info,id)
{
        //r_target = 0; seed = 2;
        generator = std::mt19937(seed);
        uniformDistribution = std::uniform_real_distribution<double>(0.0,1.0);
}

void PoissonNeuronPop::advect(std::vector<double> * synaptic_dV)
{
    // double dt           = info->dt;

    ClearSpiker();
    if (inputDependant){
    //#pragma omp parallel for
        for(unsigned long i = 0 ; i < noNeurons; i++){
            // set target rate
            lambda = synaptic_dV->at(i);
            //Check if neuron fires
            if (uniformDistribution(generator) < lambda){
                spiker.push_back(i);
            }
        }
    } else {
        int totalFiringNeurons{ binomialDistribution(generator) };
        spiker.resize(totalFiringNeurons);
        std::sample(neuronIds.begin(), neuronIds.end(), spiker.begin(),totalFiringNeurons,generator);
    }
}

void PoissonNeuronPop::LoadParameters(std::vector<std::string> *input)
{
    NeuronPop::LoadParameters(input);

    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);
        if(name.find("r_target") != std::string::npos) {
            r_target = std::stod(values.at(0));
            inputDependant = false;
            lambda 	= r_target*this->info->dt;
        } else if (name.find("seedPoisson") != std::string::npos) {
            seed = std::stoi(values.at(0));
        }
    }

    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        seed = distribution(info->globalGenerator);
        generator = std::mt19937(seed);
    }
    neuronIds.resize(noNeurons);
    std::iota(neuronIds.begin(), neuronIds.end(), 0);
    binomialDistribution=std::binomial_distribution<>(noNeurons, lambda);
}


void PoissonNeuronPop::SaveParameters(std::ofstream * stream)
{

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
    *stream <<  "#\t\tPoisson neuron: produces Poisson spiking with rate r_target (defined under stimulus). ZERO DOES NOT REMOVE THE FEATURE, YOU MUST REMOVE THE ENTIRE LINE \n";
    *stream <<  "#\t\tIf r_target not set in parameters, the neurons will fire with probability equal to the membrane potential. If Vm > 1mV, p=1 \n";

}
