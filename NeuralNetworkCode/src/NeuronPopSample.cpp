//
//  NeuronPopSample.cpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 20.11.17.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#include "NeuronPopSample.hpp"
#include "NeuronPop/HeterosynapticNeuronPop/HeteroLIFNeuronPop.hpp"


NeuronPopSample::NeuronPopSample(std::vector<std::string> *input,GlobalSimInfo * info){

    this->generalNeuronSeed = -1;
    this->info = info;
    noPopulations = 0;

    LoadParameters(input);

}

void NeuronPopSample::LoadParameters(std::vector<std::string> *input){

    std::string              name,token,type,filterStr;
    std::vector<std::string> values,neurons_strs;
    bool                     newFormat = false;
    std::vector<std::string> Ni;

    //for (std::vector<std::string>::const_iterator i = (*input).begin(); i != (*input).end(); ++i)
    //std::cout << *i << " \n";

    //Get number of populations
    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);

        if(name.find("noPopulations") != std::string::npos){
            noPopulations = static_cast<int>(std::stod(values.at(0)));
            newFormat     = true;
        }
        else if(name.find("Ni") != std::string::npos){//Deprecated. I do not even know what this is
            Ni            = values;
            noPopulations = static_cast<int>(values.size());
        }
        else if(name.find("generalNeuronSeed") != std::string::npos){
            generalNeuronSeed = std::stoi(values.at(0));
        }
    }
    neuronPops = new NeuronPop*[noPopulations];

    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        generalNeuronSeed = distribution(info->globalGenerator);
    }

    //iterate through populations
    for(unsigned int p = 0; p < noPopulations;p++){
        //std::cout << "number of population: " << std::to_string(p) << "\n";

        //filter input
        if(newFormat)
            filterStr = "neurons_" + std::to_string(p);
        else
            filterStr = "neurons_";

        FilterStringVector(input, filterStr, &neurons_strs);

        //find type
        for(std::vector<std::string>::iterator it = (*(&neurons_strs)).begin(); it != (*(&neurons_strs)).end(); ++it) {

            SplitString(&(*it),&name,&values);
            if(name.find("type") != std::string::npos){
                //std::cout << values.at(0) << "\n";
                type = values.at(0);
                break;
            }
        }

        //std::cout << "neuron type: " << type << "\n";

        //define neuronPopulation according to type
        if(type == str_LIFNeuron){
            //delete statements have been removed
            neuronPops[p] = new LIFNeuronPop(info,p);
        } else if (type == str_EIFNeuron) {
			neuronPops[p] = new EIFNeuronPop(info, p);
        } else if (type == str_PoissonNeuron) {
            neuronPops[p] = new PoissonNeuronPop(info, p);
        } else if (type == str_QIFNeuron){
            neuronPops[p] = new QIFNeuronPop(info, p);
        } else if (type == str_HeteroLIFNeuron) {
            // using more than 1 population wouSld mean there must be communcation between Connectivity objects that share the same target popluation
            neuronPops[p] = new HeteroLIFNeuronPop(info, p);
        } else if (type == str_HeteroPoissonNeuron) {
            // using more than 1 population wouSld mean there must be communcation between Connectivity objects that share the same target popluation
            neuronPops[p] = new HeteroPoissonNeuronPop(info, p);
        }else if (type == str_DictatNeuron) {
            // using more than 1 population wouSld mean there must be communcation between Connectivity objects that share the same target popluation
            neuronPops[p] = new DictatNeuronPop(info, p);
        }else {
            throw std::runtime_error(">>>Undefined type of NeuronPop.\n>>>Check Parameters.txt.");
        }

        //load parameters
        if(newFormat)
            neuronPops[p]->LoadParameters(&neurons_strs);
        else
            neuronPops[p]->LoadParameters(&neurons_strs, static_cast<unsigned long>(std::stod(Ni.at(p))));
    }

    //Set seeds
    if(generalNeuronSeed >= 0){
        std::mt19937 generator;
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        generator = std::mt19937(generalNeuronSeed);

        for(unsigned int p = 0; p < noPopulations; p++)
            neuronPops[p]->SetSeeds(distribution(generator),distribution(generator));
    }


    info->N = 0;
    for(unsigned int i = 0;i<noPopulations;i++) {
        info->N  += neuronPops[i]->GetNoNeurons();
    }

	if (info->Dimensions == 1) {
		info->Lx = info->N / static_cast<double>(info->density);
		info->Ly = 0;
	}
	else if (info->Dimensions == 2) {
		info->Lx = sqrt(info->N / static_cast<double>(info->density) );
		info->Ly = info->Lx;
	}

	for (unsigned int p = 0;p < noPopulations;p++)
		neuronPops[p]->SetPosition(neuronPops[p]->GetNoNeurons());

    // Set global neuron id
    //long global_neuron_id = 0;
    //for(int p = 0; p < noPopulations; p++){
    //    neuronPops[p]->SetGlobalNeuronId(global_neuron_id);
    //   global_neuron_id += neuronPops[p]->GetNoNeurons();
    //}
}

void NeuronPopSample::advect(std::vector<std::vector<double>> * synaptic_dV){


    for(unsigned int p = 0; p < noPopulations; p++)
        neuronPops[p]->advect(&synaptic_dV->at(p));
}


void NeuronPopSample::SaveParameters(std::ofstream * stream){

    *stream <<  "#***********************************************\n";
    *stream <<  "#************** Neuron Parameters **************\n";
    *stream <<  "#***********************************************\n";

    *stream <<  "neurons_noPopulations                 " << std::to_string(noPopulations) << "\n";
    if(info->globalSeed == -1){
        *stream <<  "neurons_generalNeuronSeed             " << std::to_string(generalNeuronSeed) << "\n";
        *stream <<  "#generalNeuronSeed = -1: seeds are defined at individual population level.\n";
        *stream <<  "#generalNeuronSeed >= 0: general seed overrides individual seeds.\n";
    }

    for(unsigned int p = 0;p < noPopulations; p++)
        neuronPops[p]->SaveParameters(stream);
}
