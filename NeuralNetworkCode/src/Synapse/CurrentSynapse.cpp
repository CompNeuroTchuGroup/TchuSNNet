#include "CurrentSynapse.hpp"

CurrentSynapse::CurrentSynapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo* infoGlobal):Synapse(targetPop,sourcePop,infoGlobal){

}

std::vector<double> CurrentSynapse::AdvectSpikers (NeuronInt spiker) {
    NeuronInt noTargets{GetNoTargetedNeurons(spiker)};
    std::vector<double> currents(noTargets,0.0);
    for(NeuronInt targetNeuronIndex : std::ranges::views::iota (0, noTargets)){
        double calcCurrent =  GetCouplingStrength(targetNeuronIndex, spiker); //If confused about syntax, talk to Antoni
        currents.at(targetNeuronIndex)+=calcCurrent;
        this->cumulatedDV   += calcCurrent;
    }
    return currents;
}

std::string CurrentSynapse::GetDataHeader(int dataColumn) {
    return "#" + std::to_string(dataColumn) + " J_"+ GetIdStr() + " (mV)\n";
}

std::string CurrentSynapse::GetUnhashedDataHeader() const {
	return "J_" + GetIdStr() + "\t";
}

std::vector<double> CurrentSynapse::GetSynapticState(NeuronInt sourceNeuron) const {
    std::vector<double> value(1);
    double Jsum = 0;
    // get average coupling strength
    for(NeuronInt targetNeuron : std::ranges::views::iota(0,GetNoTargetedNeurons(sourceNeuron))){
        Jsum += GetDistributionJ(targetNeuron,sourceNeuron);
    }
    value.at(0) = Jsum/static_cast<double>(this->GetNoTargetedNeurons(sourceNeuron));
    //value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
    return value;
}


void CurrentSynapse::LoadParameters(const std::vector<FileEntry>& synapseParameters){
    Synapse::LoadParameters(synapseParameters);
}


void CurrentSynapse::SaveParameters(std::ofstream& wParameterStream,std::string idString) const{
    Synapse::SaveParameters(wParameterStream,idString);
}
