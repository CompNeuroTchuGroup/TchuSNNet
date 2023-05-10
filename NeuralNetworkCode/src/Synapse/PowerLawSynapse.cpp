#include "PowerLawSynapse.hpp"


PowerLawSynapse::PowerLawSynapse(NeuronPop * postNeurons, NeuronPop * preNeurons, GlobalSimInfo * info) :Synapse(postNeurons, preNeurons, info)
{
	n = 1;
	Naveraging = 2;

	spike_count.resize(GetNoNeuronsPre());
	ISI_table.resize(GetNoNeuronsPre());
}


void PowerLawSynapse::advectSpikers(std::vector<double>& currents, long spiker)
{
	double ISI;
	double nu_H; //estimated firing rate based on ISI
	// double t = info->time_step*info->dt;
	int N;

	spike_count[spiker] = spike_count[spiker]+1;
	ISI_table[spiker][(spike_count[spiker] - 1) % Naveraging] = neuronsPre->GetTimeSinceLastSpike(spiker);
	N = std::min(spike_count[spiker], Naveraging);
	ISI = ISI_table[spiker].sum()/N;//Average ISI over the last Naveraging spikes
	nu_H = 1 / ISI;
	if (N > 1)
		nu_H = static_cast<double>(N - 1) / N * nu_H;//the expected value of 1/ISI is nu*N/(N-1) for a Poisson process
	std::vector<unsigned long> *tL = geometry->GetTargetList(spiker);

	for (unsigned int target = 0; target < tL->size(); target++) {
		double c = GetCouplingStrength(spiker, target);
		currents[target] += c* pow(nu_H, n);
		this->cumulatedDV += c * pow(nu_H, n);
	}
}


void PowerLawSynapse::LoadParameters(std::vector<std::string> *input) {

	Synapse::LoadParameters(input);

	std::string              name;
	std::vector<std::string> values;

	for (auto & it : *input) {
		SplitString(&it, &name, &values);

		if (name.find("powerlawsyn_n") != std::string::npos) {
			n = std::stod(values.at(0));
		}
		if (name.find("powerlawsyn_N_averaging") != std::string::npos) {
			Naveraging = std::stoi(values.at(0));
		}
	}
	for (unsigned long source = 0;source < GetNoNeuronsPre();source++) {
		ISI_table[source].resize(Naveraging);
		spike_count[source] = 0;
	}

}

void PowerLawSynapse::SaveParameters(std::ofstream * stream, std::string id_str) {
	Synapse::SaveParameters(stream, id_str);
	*stream << id_str << "powerlawsyn_n\t\t\t\t\t" << std::to_string(n) << "\n";
	*stream << id_str << "powerlawsyn_N_averaging\t\t" << std::to_string(Naveraging) << "\n";
	*stream << "#\t\tSynaptic strength as a power law of the presynaptic firing rate : (J_eff=J/<ISI>^n); with the mean of the last N_averaging ISIs measured in seconds\n";
}


std::string PowerLawSynapse::GetDataHeader(int data_column)
{
	return "#" + std::to_string(data_column) + " J" + GetIdStr() + " Synaptic weight at 1Hz \n";
}

std::string PowerLawSynapse::GetUnhashedDataHeader() {
	return "J_" + GetIdStr() + "\t";
}

std::valarray<double> PowerLawSynapse::GetSynapticState(int pre_neuron){
	std::valarray<double> val(1);
	double Jsum = 0;
	// get average coupling strength
	for (unsigned int target = 0; target < this->GetNumberOfPostsynapticTargets(pre_neuron); target++) {
		Jsum += *(geometry->GetDistributionJ(pre_neuron, target));
	}
	val[0] = Jsum / double(this->GetNumberOfPostsynapticTargets(pre_neuron));
	//val[0] = GetCouplingStrength()*double(this->GetNumberOfPostsynapticTargets(pre_neuron));
	return val;
}
