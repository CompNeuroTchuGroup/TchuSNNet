#include "PModelSynapse.hpp"

PModelSynapse::PModelSynapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal) : Synapse(targetPop, sourcePop, infoGlobal) {
    hasPlasticityModel=true;
}
std::vector<double> PModelSynapse::AdvectSpikers(NeuronInt spiker) {
    NeuronInt noTargets{GetNoTargetedSynapses(spiker)};
    std::vector<double> currents(noTargets,0.0);
    for(NeuronInt targetNeuronIndex : std::ranges::views::iota (0, noTargets)){
    // for(NeuronInt targetNeuronIndex{0};targetNeuronIndex<noTargets;targetNeuronIndex++){
        this->targetPop->RecordExcitatorySynapticSpike(targetSpineList.at(spiker).at(targetNeuronIndex).first, targetSpineList.at(spiker).at(targetNeuronIndex).second->idInMorpho);
        currents.at(targetNeuronIndex)=targetSpineList.at(spiker).at(targetNeuronIndex).second->GetWeight();
    }
    return currents;
}

void PModelSynapse::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
    Synapse::SaveParameters(wParameterStream, idString);
    
    targetPop->SavePlasticityModel(wParameterStream, idString+"pmodel_");//!!!!

    targetPop->ModelSaved();
    if (geometry != nullptr){
        wParameterStream << "########### Connectivity ################\n";
        geometry->SaveParameters(wParameterStream,idString);
    }
}

void PModelSynapse::LoadParameters(const std::vector<FileEntry> &hybridParameters) {
    targetPop->LoadPlasticityModel(FilterStringEntries(hybridParameters, "pmodel"));
    this->ignoreJDParameters=targetPop->ignoreJDParameters();
    Synapse::LoadParameters(hybridParameters);
    Dmax=targetPop->GetMaxGapDelay(delayPerMicrometer);
}
void PModelSynapse::AllocateSynapse(NeuronInt targetNeuron, NeuronInt sourceNeuron) {
    AllocateSynapseWithPlasticity(targetNeuron, sourceNeuron);
}
std::string PModelSynapse::GetDataHeader(int dataColumn) {
    return "#" + std::to_string(dataColumn) + " J_"+ GetIdStr() + " (mV)\n";
}

std::string PModelSynapse::GetUnhashedDataHeader() const {
	return "J_" + GetIdStr() + "\t";
}

std::vector<double> PModelSynapse::GetSynapticState(NeuronInt sourceNeuron) const {
    std::vector<double> value(1);
    double Jsum = 0;
    // get average coupling strength
    Jsum = std::accumulate(targetSpineList.at(sourceNeuron).begin(), targetSpineList.at(sourceNeuron).end(), 0.0, [](double Jsum, std::pair<NeuronInt, BaseSpinePtr> synapse){
        return Jsum + synapse.second->GetWeight();
    });
    value.at(0) = Jsum/static_cast<double>(this->GetNoTargetedSynapses(sourceNeuron));
    //value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
    return value;
}

