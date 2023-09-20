#include "PModelSynapse.hpp"

PModelSynapse::PModelSynapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal) : Synapse(targetPop, sourcePop, infoGlobal) {
    hasPlasticityModel=true;
}
std::vector<double> PModelSynapse::AdvectSpikers(NeuronInt spiker) {
    NeuronInt noTargets{GetNoTargetedNeurons(spiker)};
    std::vector<double> currents(noTargets,0.0);
    for(NeuronInt targetNeuronIndex : std::ranges::views::iota (0, noTargets)){
    // for(NeuronInt targetNeuronIndex{0};targetNeuronIndex<noTargets;targetNeuronIndex++){
        double current =  targetSpineList.at(spiker).at(targetNeuronIndex).second->GetWeight(); //If confused about syntax, talk to Antoni

        this->targetPop->RecordExcitatorySynapticSpike(targetSpineList.at(spiker).at(targetNeuronIndex).first, targetSpineList.at(spiker).at(targetNeuronIndex).second->GetIdInMorpho());
        
        currents.at(targetNeuronIndex)+=current;
        this->cumulatedDV   += current;
    }
    return currents;
}

void PModelSynapse::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
    wParameterStream << idString << "type\t\t\t\t\t\t" << GetTypeStr() << "\n";
    wParameterStream << idString << "connected\t\t\t\t\t\t" << std::boolalpha << this->isConnectedBool <<std::noboolalpha<< "\n";
    wParameterStream << idString << "D_min\t\t\t\t\t\t" << std::to_string(this->Dmin*infoGlobal->dtTimestep) << " #secs\n";
    wParameterStream << idString << "D_max\t\t\t\t\t\t" << std::to_string(this->Dmax*infoGlobal->dtTimestep) << " #secs\n";
    if (!this->ignoreJDistribution){
        wParameterStream << idString << "J\t\t\t\t\t\t\t" << std::to_string(this->J) << " #dmV/Spike\n";
        wParameterStream << idString << "Sigma_j\t\t\t\t\t\t" << std::to_string(this->sigmaJ) << " #dmV/Spike\n";
        wParameterStream << idString << "J_pot\t\t\t\t\t\t" << std::to_string(this->Jpot) << " #dmV/Spike\n";
        wParameterStream << idString << "P_pot\t\t\t\t\t\t" << std::to_string(this->Ppot) << "\n";
    }
    if (userSeed){
        wParameterStream << idString << "seed\t\t\t\t\t\t" << std::to_string(this->seed)  << "\n";
    }

    //Synaptic targeting
    if (this->targetPop->IsBranchedBool()){
        wParameterStream << idString << "targetBranch\t\t\t\t\t";
        if (this->branchTarget.randomTargetBranch){
            wParameterStream<<"random"; //Missing comments on what this is supposed to do
        } else if(this->branchTarget.setTargetBranch){
            wParameterStream<<std::to_string(this->branchTarget.targetBranch);//Missing comments on what this is supposed to do
        }else if(this->branchTarget.orderedTargetBranch){
        wParameterStream<<"ordered";
        } else {
            wParameterStream<<"none";//Missing comments on what this is supposed to do
        }
        wParameterStream << "\t\t\t#You can target branches in an 'ordered' manner (0,1,2...), 'random', or set (if you input a number). Put none if the HS does not used branched morphology\n";
        wParameterStream << idString+"pmodel_"<<"slotOrder\t\t";
        if (this->branchTarget.firstSlotTrueLastSlotFalse){
            wParameterStream << "first\t";
        } else {
            wParameterStream<<"last\t";
        }
        wParameterStream << "\t"<<"#'first' synapse allocation will allocate synapses from the beggining to the end of the available slots. 'last' will do the opposite. This only makes sense in ordered allocation\n";

    }
    if (this->ignoreJDistribution){
        wParameterStream << idString << "relativeCoupling\t\t\t\t\t" << std::to_string(relativeCouplingStrength)<< "\t#Relative coupling strength (only non-J plasticity models)\n";
    }
    targetPop->SavePlasticityModel(wParameterStream, idString+"pmodel_");//!!!!

    targetPop->ModelSaved();
    if (geometry != nullptr){
        wParameterStream << "########### Connectivity ################\n";
        geometry->SaveParameters(wParameterStream,idString);
    }
}

void PModelSynapse::LoadParameters(const std::vector<FileEntry> &hybridParameters) {
    targetPop->LoadPlasticityModel(FilterStringEntries(hybridParameters, "pmodel"));
    this->ignoreJDistribution=targetPop->ignoreJDistribution();
    Synapse::LoadParameters(hybridParameters);
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
    value.at(0) = Jsum/static_cast<double>(this->GetNoTargetedNeurons(sourceNeuron));
    //value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
    return value;
}

