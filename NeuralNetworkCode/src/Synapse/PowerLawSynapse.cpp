#include "PowerLawSynapse.hpp"


PowerLawSynapse::PowerLawSynapse(PopPtr targetPop, PopPtr sourceNeurons, GlobalSimInfo* infoGlobal) :Synapse(targetPop, sourceNeurons, infoGlobal) {
	spikeCountVector.resize(GetNoSourceNeurons(), 0);
	tableISI.resize(GetNoSourceNeurons());
}


std::vector<double> PowerLawSynapse::AdvectSpikers (NeuronInt spiker) {
	double ISI;
	double averageFrequency; //estimated firing rate based on ISI (previously nu_H)
	// double t = infoGlobal->timeStep*infoGlobal->dtTimestep;
	int noSpikes;
	std::vector<double>& singleISITable = tableISI.at(spiker);

	singleISITable.at((spikeCountVector.at(spiker)) % noSpikesForFiringRate) = sourcePop->GetTimeSinceLastSpike(spiker);
	spikeCountVector.at(spiker) += 1;
	
	noSpikes = std::min(spikeCountVector.at(spiker), noSpikesForFiringRate);
	// ISI = tableISI.at(spiker).sum()/noSpikes;//Average ISI over the last Naveraging spikes
	ISI = std::reduce(singleISITable.begin(), singleISITable.end()) /noSpikes;//Average ISI over the last Naveraging spikes. Reduce does the .sum()
	averageFrequency = 1 / ISI;
	if (noSpikes > 1){
		averageFrequency = (static_cast<double>(noSpikes - 1) / noSpikes) * averageFrequency;//the expected value of 1/ISI is nu*N/(N-1) for a Poisson process
	}
	double exponent{pow(averageFrequency, kExponent_PLaw)};
	std::vector<double> currents(targetSpineList.at(spiker).size());
	for (NeuronInt targetNeuron : std::ranges::views::iota(0, static_cast<NeuronInt>(targetSpineList.at(spiker).size()))) {
		double couplingStrength = GetCouplingStrength(targetNeuron, spiker)* exponent;
		currents.at(targetNeuron) = couplingStrength;//No plus equal because the currents vector is based on indexes, not actual neuron ids
		this->cumulatedDV += couplingStrength;
	}
	return currents;
}


void PowerLawSynapse::LoadParameters(const std::vector<FileEntry>& synapseParameters) {

	Synapse::LoadParameters(synapseParameters);

	for (auto& [parameterName, parameterValues] : synapseParameters) {
		if (parameterName.find("powerlawsyn_n") != std::string::npos) {
			kExponent_PLaw = std::stod(parameterValues.at(0));
		} else if (parameterName.find("powerlawsyn_N_averaging") != std::string::npos) {
			noSpikesForFiringRate = std::stoi(parameterValues.at(0));
		}
	}
	for (NeuronInt sourceNeuron : std::ranges::views::iota(0,GetNoSourceNeurons())) {
		tableISI.at(sourceNeuron).resize(noSpikesForFiringRate);
		// spikeCountVector.at(sourceNeuron) = 0;
	}

}

void PowerLawSynapse::SaveParameters(std::ofstream& wParameterStream, std::string idString) const {
	Synapse::SaveParameters(wParameterStream, idString);
	wParameterStream << idString << "powerlawsyn_n\t\t\t\t\t" << std::to_string(kExponent_PLaw) << "\n";
	wParameterStream << idString << "powerlawsyn_N_averaging\t\t\t\t" << std::to_string(noSpikesForFiringRate) << "\n";
	wParameterStream << "#\t\tSynaptic strength as a power law of the presynaptic firing rate : (J_eff=J/<ISI>^n); with the mean of the last N_averaging ISIs measured in seconds\n";
}


std::string PowerLawSynapse::GetDataHeader(int data_column) {
	return "#" + std::to_string(data_column) + " J" + GetIdStr() + " Synaptic weight at 1Hz \n";
}

std::string PowerLawSynapse::GetUnhashedDataHeader() const {
	return "J_" + GetIdStr() + "\t";
}

std::vector<double> PowerLawSynapse::GetSynapticState(NeuronInt sourceNeuron) const{
	std::vector<double> synapticState(1);
	double Jsum = 0;
	// get average coupling strength
	for (NeuronInt targetNeuron : std::ranges::views::iota(0, this->GetNoTargetedNeurons(sourceNeuron))) {
		Jsum += GetDistributionJ(targetNeuron,sourceNeuron);
	}
	synapticState.at(0) = Jsum / static_cast<double>(this->GetNoTargetedNeurons(sourceNeuron));
	//value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
	return synapticState;
}
