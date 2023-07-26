
#include "Stimulus.hpp"

Stimulus::Stimulus(std::shared_ptr<NeuronPopSample>  neuronSample,GlobalSimInfo*  infoGlobal) : infoGlobal{infoGlobal}, neurons{neuronSample} {
    //int N         = neurons->GetTotalNeurons();
    //signal_array = new double*[totalNeuronPops]//Or sth like this
    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
        signalMatrix.push_back(std::vector<double>(neurons->GetNeuronsPop(neuronPop), 0.0));
    }

    // for(PopInt neuronPop = 0; neuronPop < totalNeuronPops; neuronPop++){
    //     signal_array.at(neuronPop) = new double[neurons->GetNeuronsPop(neuronPop)];
    //     for(signed long i = 0; i < neurons->GetNeuronsPop(totalNeuronPops); i++)
    //         signal_array[totalNeuronPops][i] = 0.0;
    // }
    std::uniform_int_distribution<int> distribution(0,INT_MAX);
    generator=std::mt19937(distribution(infoGlobal->globalGenerator));
}

void Stimulus::SaveParameters(std::ofstream& wParameterStream) const{
    wParameterStream <<  "#**************************************************\n";
    wParameterStream <<  "#************** Stimulus Parameters ***************\n";
    wParameterStream <<  "#**************************************************\n";
    wParameterStream <<  "stimulus_type                        " << GetType() << "\n";
    if (userSeed){
        wParameterStream << "stimulus_seed\t\t\t\t" << std::to_string(this->seed)  << "\n";
    }

}

void Stimulus::LoadParameters(const std::vector<FileEntry>& stimulusParameters) {
    for(auto& [parameterName, parameterValues] : stimulusParameters){
        if(parameterName.find("seed") != std::string::npos){
            seed=std::stoi(parameterValues.at(0));
            generator = std::mt19937(seed);
            userSeed=true;
        }
    }
}

void Stimulus::Update(std::vector<std::vector<double>>& synaptic_dV){
    for(PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
		for (NeuronInt neuron : std::ranges::views::iota(0, neurons->GetNeuronsPop(neuronPop)))
			synaptic_dV.at(neuronPop).at(neuron) += signalMatrix.at(neuronPop).at(neuron);
    }
}
