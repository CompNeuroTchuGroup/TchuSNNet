#include "PoissonNeuronPop.hpp"

PoissonNeuronPop::PoissonNeuronPop(GlobalSimInfo* infoGlobal,NeuronInt neuronID):NeuronPop(infoGlobal,neuronID) {
        //targetRate = 0; seed = 2;
        // generator = std::mt19937(seed);
        uniformDistribution = std::uniform_real_distribution<double>(0.0,1.0);
}

void PoissonNeuronPop::Advect(const std::vector<double>& synaptic_dV) {
    // double dtTimestep           = infoGlobal->dtTimestep;
    ClearSpikerVector();
    if (inputDependant){
    //#pragma omp parallel for
        for(NeuronInt neuron : std::ranges::views::iota(0, noNeurons)){
            // set target rate
            lambda = synaptic_dV.at(neuron);
            //Check if neuron fires
            if (uniformDistribution(generator) < lambda){//this is essentially bernoulli trial, we call Poisson because low prob of 1, negligible of >1 and small timesteps
                spikerNeurons.push_back(neuron);
            }
        }
    } else {
        NeuronInt totalFiringNeurons{binomialDistribution(generator)};
        spikerNeurons.resize(totalFiringNeurons);
        std::sample(neuronIds.begin(), neuronIds.end(), spikerNeurons.begin(),totalFiringNeurons,generator);
    }
    this->AdvectPlasticityModel();
}

void PoissonNeuronPop::LoadParameters(const std::vector<FileEntry>& neuronParameters) {
    NeuronPop::LoadParameters(neuronParameters);

    for(auto& [parameterName, parameterValues] : neuronParameters) {
        if(parameterName.find("r_target") != std::string::npos) {
            targetRate = std::stod(parameterValues.at(0));
            inputDependant = false;
            lambda 	= targetRate*this->infoGlobal->dtTimestep;
        // } else if (parameterName.find("seedPoisson") != std::string::npos) {
        //     seed = std::stoi(parameterValues.at(0));
        }
    }
    if (!inputDependant && hasPlasticity){
        throw "PoissonNeuron may not have a plasticity model if it is not input dependant";
    }

    // if(infoGlobal->globalSeed != -1){
    //     std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    //     seed = distribution(infoGlobal->globalGenerator);
    //     generator = std::mt19937(seed);
    // }
    neuronIds.resize(noNeurons);
    std::iota(neuronIds.begin(), neuronIds.end(), 0);
    binomialDistribution=std::binomial_distribution<>(noNeurons, lambda);
}


void PoissonNeuronPop::SaveParameters(std::ofstream& wParameterStream) const {

    std::string idString = "neurons_" + std::to_string(GetId());

    //std::cout<<"printing neuron id: "<<id<<"\n";

    //NeuronPop::SaveParameters(stream);
    wParameterStream <<  "#***********************************************\n";
    wParameterStream <<  idString + "_noNeurons\t\t\t" << noNeurons << "\n";
    wParameterStream <<  idString + "_type\t\t\t" << GetType() << "\n";
	// if (infoGlobal->globalSeed == -1) {
	// 	*stream << id + "_seedPoisson                 " << std::to_string(seed) << "\n";
	// }
    if (!inputDependant){
        wParameterStream <<  idString + "_r_target\t\t\t" << std::to_string(targetRate)  << "\n";
    }
    if(userSeeds){
        // wParameterStream <<  neuronID + "_seedInitialPotentials   " << this->seedInitialPotentials << "\n";
        // wParameterStream <<  neuronID + "_seedInitialPrevSpike    " << this->seedInitialPreviousSpike << "\n";
        wParameterStream <<  idString + "_seed\t\t\t\t" << this->seed << "\n";
    }
    wParameterStream <<  "#\t\tPoisson neuron: produces Poisson spiking with rate r_target (defined under stimulus). ZERO DOES NOT REMOVE THE FEATURE, YOU MUST REMOVE THE ENTIRE LINE \n";
    wParameterStream <<  "#\t\tIf r_target not set in parameters, the neurons will fire with probability equal to the membrane potential. If Vm > 1mV, p=1 \n";

}
