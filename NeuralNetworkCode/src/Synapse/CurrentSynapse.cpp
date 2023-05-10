#include "CurrentSynapse.hpp"

CurrentSynapse::CurrentSynapse(NeuronPop * postNeurons,NeuronPop * preNeurons,GlobalSimInfo * info):Synapse(postNeurons,preNeurons,info) {}

void CurrentSynapse::advectSpikers (std::vector<double>& currents, long spiker)
{
    std::vector<unsigned long> *tL = geometry->GetTargetList(spiker);

    for(unsigned int target = 0; target < tL->size(); target++){
        double c = GetCouplingStrength(spiker, target);
        currents[target] += c;
        this->cumulatedDV   += c;
    }

}

std::string CurrentSynapse::GetDataHeader(int data_column)
{
    return "#" + std::to_string(data_column) + " J_"+ GetIdStr() + " (mV)\n";
}

std::string CurrentSynapse::GetUnhashedDataHeader()
{
	return "J_" + GetIdStr() + "\t";
}

std::valarray<double> CurrentSynapse::GetSynapticState(int pre_neuron)
{
    std::valarray<double> val(1);
    double Jsum = 0;
    // get average coupling strength
    for(unsigned int target=0; target < this->GetNumberOfPostsynapticTargets(pre_neuron); target++){
        Jsum += *(geometry->GetDistributionJ(pre_neuron,target));
    }
    val[0] = Jsum/double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    //val[0] = GetCouplingStrength()*double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    return val;
}


void CurrentSynapse::LoadParameters(std::vector<std::string> *input){
    Synapse::LoadParameters(input);
}


void CurrentSynapse::SaveParameters(std::ofstream * stream,std::string id_str){
    Synapse::SaveParameters(stream,id_str);
}
