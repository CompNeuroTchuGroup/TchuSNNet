#include "HeteroLIFNeuronPop.hpp"

HeteroLIFNeuronPop::HeteroLIFNeuronPop(GlobalSimInfo * info, int id): HeteroNeuronPop(info, id) {
    reset_type = 0;
}

void HeteroLIFNeuronPop::advect(std::vector<double> * synaptic_dV)
{
    ClearSpiker();

    //#pragma omp parallel for
    for(unsigned long i = 0 ; i < noNeurons; i++)
    {
        //(a) wait for refractory period
        if((info->time_step-previous_spike_step[i]) <= refractorySteps)
            continue;

        //(b) advect
        potential[i] = potential[i]*membraneExpDecay + synaptic_dV->at(i);

        //(c) determine if neuron has spiked
        if(potential[i] > v_thresh)
        {
            spiker.push_back(i);

            //Reset potential
            if(reset_type == 0)
                potential[i] = v_reset;
            else if(reset_type == 1){
                while(potential[i] > v_thresh)
                    potential[i] = v_reset + (potential[i] - v_thresh);
            }

        }
    }

    for (auto neuron: this->spiker) {
        this->morphology[neuron]->RecordPostSpike();
    }

    for (unsigned long morphId =  0; morphId < morphology.size(); ++morphId) {
        this->morphology[morphId]->advect();
    }

}

void HeteroLIFNeuronPop::LoadParameters(std::vector<std::string> *input){

    HeteroNeuronPop::LoadParameters(input);

    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);
        if(name.find("resetType") != std::string::npos){
            reset_type = std::stoi(values.at(0));
        }
    }
    membraneExpDecay=exp(-info->dt/tau_m);
}

void HeteroLIFNeuronPop::SaveParameters(std::ofstream * stream){

    std::string id = "neurons_" + std::to_string(GetId());

    HeteroNeuronPop::SaveParameters(stream);
    *stream <<  id + "_resetType                   " << std::to_string(reset_type)  << "\n";
    *stream <<  "#\t\tHeteroLIF neuron: dV/dt = -V/tau_m + RI/tau_m \n";
    *stream <<  "#\t\tresetType 0: v = v_reset\n";
    *stream <<  "#\t\tresetType 1: v = v_reset + (v - v_thresh) \n";

}

std::string HeteroLIFNeuronPop::GetType() {
    return str_HeteroLIFNeuron;
}
