//
//  NeuronPopSample.cpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 20.11.17.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#include "NeuronPopSample.hpp"


NeuronPopSample::NeuronPopSample(std::vector<FileEntry> neuronParameters, GlobalSimInfo* infoGlobal) : infoGlobal{infoGlobal}{

    LoadParameters(neuronParameters);
}

void NeuronPopSample::LoadParameters(const std::vector<FileEntry>& neuronParameters){

    //Get number of populations
    for(auto& [parameterName, parameterValues] : neuronParameters) {

        if(parameterName.find("noPopulations") != std::string::npos){
            noPopulations = static_cast<PopInt>(std::stoi(parameterValues.at(0)));
        }
        // else if(parameterName.find("generalNeuronSeed") != std::string::npos){
        //     generalNeuronSeed = std::stoi(parameterValues.at(0));
        // }
    }
    //neuronPops.resize(noPopulations);

    // if(infoGlobal->globalSeed != -1){
    //     std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    //     generalNeuronSeed = distribution(infoGlobal->globalGenerator);
    // }

    //iterate through populations
    for(PopInt population : std::ranges::views::iota(0,noPopulations)){
        //It could be argued that the previous implementation was somewhat safer to misallocation of neuronPop ids, but the single-point of failure (the population signed int in the for loop statement),
        //is exactly the same. The id that shows in neurons_id is the end-all-be-all of the position (which follows ascending order always).
        //std::cout << "number of population: " << std::to_string(p) << "\n";
        std::string neuronType;
        //filter population parameters
        std::string filterStr = "neurons_" + std::to_string(population);

        std::vector<FileEntry> singleNeuronPopParameters = FilterStringEntries(neuronParameters, filterStr);

        //find type
        for(FileEntry singleEntry : singleNeuronPopParameters) {
            if (singleEntry.parameterNameContains("type") && !singleEntry.parameterNameContains("morphology")){
                LoadNeuronPop(singleEntry.parameterValues.at(0), population, singleNeuronPopParameters);
                break;
            }
        }
    }
    std::for_each(neuronPops.begin(), neuronPops.end(), [this](const PopPtr& neuronPop){
        infoGlobal->totalNeurons  += neuronPop->GetNoNeurons();
    });
    // for(PopPtr neuronPopPtr : neuronPops) {
    //     infoGlobal->totalNeurons  += neuronPop->GetNoNeurons();
    // }
	if (infoGlobal->dimensions == 1) {
		infoGlobal->xAxisLength = infoGlobal->totalNeurons / static_cast<double>(infoGlobal->density);
	} else if (infoGlobal->dimensions == 2) {
		infoGlobal->xAxisLength = sqrt(infoGlobal->totalNeurons / static_cast<double>(infoGlobal->density) );
		infoGlobal->yAxisLength = infoGlobal->xAxisLength;
	}
    std::for_each(neuronPops.begin(), neuronPops.end(), [](const PopPtr& neuronPop){
        neuronPop->SetPosition();
    });
	// for (PopInt neuronPop : std::ranges::views::iota(0, noPopulations)){
	// 	neuronPops.at(neuronPop)->SetPosition(neuronPops.at(neuronPop)->GetNoNeurons());
    // }

}

void NeuronPopSample::LoadNeuronPop(std::string neuronPopType, PopInt popID, std::vector<FileEntry> neuronParameters) {

        if(neuronPopType == IDstringLIFNeuron || neuronPopType ==IDstringHeteroLIFNeuron){
            neuronPops.push_back(std::make_shared<LIFNeuronPop>(infoGlobal,popID));
        } else if (neuronPopType == IDstringEIFNeuron) {
			neuronPops.push_back(std::make_shared<EIFNeuronPop>(infoGlobal, popID));
        } else if (neuronPopType == IDstringPoissonNeuron || neuronPopType == IDstringHeteroPoissonNeuron) {
            neuronPops.push_back(std::make_shared<PoissonNeuronPop>(infoGlobal, popID));
        } else if (neuronPopType == IDstringQIFNeuron){
            neuronPops.push_back(std::make_shared<QIFNeuronPop>(infoGlobal, popID));
        }else if (neuronPopType == IDstringDictatNeuron) {
            neuronPops.push_back(std::make_shared<DictatNeuronPop>(infoGlobal, popID));
        }else if (neuronPopType == IDstringCorrelatedPoissonNeuron) {
            neuronPops.push_back(std::make_shared<CorrPoissonNeuronPop>(infoGlobal, popID));
        }else {
            throw ">>>Undefined type of NeuronPop.\n>>>Check Parameters.txt.";
        }

        //load parameters
        neuronPops.back()->LoadParameters(neuronParameters);
}

void NeuronPopSample::Advect(std::vector<std::vector<double>>& synaptic_dV){
    // for(PopInt neuronPop : std::ranges::views::iota(0, noPopulations)) //Here container iteration is not possible because we need the index
    //     neuronPops.at(neuronPop)->Advect(synaptic_dV.at(neuronPop));
    // std::for_each(neuronPops.begin(), neuronPops.end(), [synaptic_dV](PopPtr neuronPopPtr){
    //     neuronPopPtr->Advect(synaptic_dV.at(neuronPopPtr->GetId()));
    // });
    std::for_each(PAR,neuronPops.begin(), neuronPops.end(), [synaptic_dV](PopPtr neuronPopPtr){
        neuronPopPtr->Advect(synaptic_dV.at(neuronPopPtr->GetId()));
    });
    //Test9 debugging
    // std::cout << "spikers:";
    // for (double spiker : neuronPops.at(0)->GetSpikers()) {
    //     std::cout << spiker;
    // }
    // std::cout << std::endl;
    /*    std::cout << "spikers:";
    int spikers{};
    for (unsigned int p = 0; p < noPopulations; p++)
        spikers += neuronPops.at(p)->GetSpikers().size();
    std::cout <<spikers<< std::endl;*/

}


void NeuronPopSample::SaveParameters(std::ofstream& wParameterStream) const {

    wParameterStream <<  "#***********************************************\n";
    wParameterStream <<  "#************** Neuron Parameters **************\n";
    wParameterStream <<  "#***********************************************\n";

    wParameterStream <<  "neurons_noPopulations                 " << std::to_string(noPopulations) << "\n";
    // if(infoGlobal->globalSeed == -1){
    //     wParameterStream <<  "neurons_generalNeuronSeed             " << std::to_string(generalNeuronSeed) << "\n";
    //     wParameterStream <<  "#generalNeuronSeed = -1: seeds are defined at individual population level.\n";
    //     wParameterStream <<  "#generalNeuronSeed >= 0: general seed overrides individual seeds.\n";
    // }

    for(PopPtr neuronPopPtr: neuronPops){
        neuronPopPtr->SaveParameters(wParameterStream);
    }
}
