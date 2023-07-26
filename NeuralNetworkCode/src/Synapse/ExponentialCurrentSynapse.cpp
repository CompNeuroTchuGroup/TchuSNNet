#include "ExponentialCurrentSynapse.hpp"

ExponentialCurrentSynapse::ExponentialCurrentSynapse(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal):Synapse(targetPop,sourcePop,infoGlobal){}

void ExponentialCurrentSynapse::ResetWaitingMatrixEntry() {
    //std::vector<double> currents;
	int currentIndex = (this->infoGlobal->timeStep) % (Dmax + 1);


	if (Dmax > 0) {
		int nextIndex = (this->infoGlobal->timeStep + 1) % (Dmax + 1);
        for (std::vector<double>& targetWaitingVector : waitingMatrix){
            targetWaitingVector.at(nextIndex) += expDecayConstant * targetWaitingVector.at(currentIndex);//exponential synapse
			targetWaitingVector.at(currentIndex) = 0;//reset
        }
		// for (signed long targetNeuron : std::ranges::views::iota(0,targetPop->GetNoNeurons())) {
		// 	waitingMatrix.at(targetNeuron).at(nextIndex) += expDecayConstant * waitingMatrix.at(targetNeuron).at(currentIndex);//exponential synapse
		// 	waitingMatrix.at(targetNeuron).at(currentIndex) = 0;//reset
		// }
	}
	else {
        for (std::vector<double>& targetWaitingVector : waitingMatrix){
            targetWaitingVector.at(currentIndex) *= expDecayConstant; //* waitingMatrix.at(targetNeuron).at(currentIndex);//reset
        }
		// for (signed long targetNeuron : std::ranges::views::iota(0,targetPop->GetNoNeurons()))
		// 	waitingMatrix.at(targetNeuron).at(currentIndex) *= expDecayConstant; //* waitingMatrix.at(targetNeuron).at(currentIndex);//reset
	}
}

void ExponentialCurrentSynapse::ResetcumulatedDV() {
	cumulatedDV *= expDecayConstant;
}

std::vector<double> ExponentialCurrentSynapse::AdvectSpikers (NeuronInt spiker){
    NeuronInt noTargets{GetNoTargetedNeurons(spiker)};
    std::vector<double> currents(noTargets);
    for(NeuronInt targetNeuron : std::ranges::views::iota (0, noTargets)){
		double calcCurrent = GetCouplingStrength(targetNeuron, spiker)*infoGlobal->dtTimestep/tauConstant;
		currents.at(targetNeuron)+=calcCurrent;
		this->cumulatedDV += calcCurrent;
	}
    return currents;
}

std::string ExponentialCurrentSynapse::GetDataHeader(int dataColumn){
    return "#" + std::to_string(dataColumn) + " J_"+ GetIdStr() +" (mV)\n";
}

std::string ExponentialCurrentSynapse::GetUnhashedDataHeader() const {
	return "J_" + GetIdStr() + "\t";
}

std::vector<double> ExponentialCurrentSynapse::GetSynapticState(NeuronInt sourceNeuron) const {
    std::vector<double> value(1);
    // get average coupling strength
    double Jsum = 0;
    for(NeuronInt targetNeuron : std::ranges::views::iota(0,this->GetNoTargetedNeurons(sourceNeuron))){
        Jsum += GetDistributionJ(targetNeuron,sourceNeuron);
    }

    value.at(0) = Jsum/static_cast<double>(this->GetNoTargetedNeurons(sourceNeuron));
    //value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
    return value;
}

void ExponentialCurrentSynapse::LoadParameters(const std::vector<FileEntry>& synapseParameters){
    Synapse::LoadParameters(synapseParameters);

    for(auto& [parameterName, parameterValues] : synapseParameters) {
        if(parameterName.find("exponential_tau") != std::string::npos){
            tauConstant=(std::stod(parameterValues.at(0)));
            expDecayConstant = exp(-infoGlobal->dtTimestep / tauConstant);
        }
    }
}


void ExponentialCurrentSynapse::SaveParameters(std::ofstream& wParameterStream,std::string idString) const{
    Synapse::SaveParameters(wParameterStream,idString);

    wParameterStream << idString << "exponential_tau\t\t\t\t\t" << std::to_string(tauConstant) << " #secs\n";
    wParameterStream << "#\t\tThe Post Synaptic Potential decays exponentially through time. The AUC is determined by J and does not depend on tau\n";
}
