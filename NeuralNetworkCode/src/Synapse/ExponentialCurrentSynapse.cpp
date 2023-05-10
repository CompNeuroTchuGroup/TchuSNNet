#include "ExponentialCurrentSynapse.hpp"

ExponentialCurrentSynapse::ExponentialCurrentSynapse(NeuronPop * postNeurons,NeuronPop * preNeurons,GlobalSimInfo * info):Synapse(postNeurons,preNeurons,info){}

//void ExponentialCurrentSynapse::advect_finalize (std::vector<double> * synaptic_dV, std::vector<std::vector<std::vector<double>>> * waiting_matrix)
void ExponentialCurrentSynapse::ResetWaitingMatrixEntry()
{
    // double		dt = info->dt;
	double		h = exp(-info->dt / Tau);
    std::vector<double> currents;
	int current_index = (this->info->time_step) % (D_max + 1);


	if (D_max > 0) {
		int next_index = (this->info->time_step + 1) % (D_max + 1);
		for (unsigned long i = 0; i < neuronsPost->GetNoNeurons(); i++) {
			waiting_matrix[i][next_index] += h * waiting_matrix[i][current_index];//exponential synapse
			waiting_matrix[i][current_index] = 0;//reset
		}
	}
	else {
		for (unsigned long i = 0; i < neuronsPost->GetNoNeurons(); i++)
			waiting_matrix[i][current_index] = h * waiting_matrix[i][current_index];//reset
	}
}

void ExponentialCurrentSynapse::resetcumulatedDV() {
	double		h = exp(-info->dt / Tau);
	cumulatedDV = cumulatedDV*h;
}

void ExponentialCurrentSynapse::advectSpikers (std::vector<double>& currents, long spiker){
	std::vector<unsigned long> *tL = geometry->GetTargetList(spiker);

	for (unsigned int target = 0; target < tL->size(); target++) {
		double c = GetCouplingStrength(spiker, target)*info->dt/Tau;
		currents[target] += c;
		this->cumulatedDV += c;
	}
}

std::string ExponentialCurrentSynapse::GetDataHeader(int data_column){
    return "#" + std::to_string(data_column) + " J_"+ GetIdStr() +" (mV)\n";
}

std::string ExponentialCurrentSynapse::GetUnhashedDataHeader() {
	return "J_" + GetIdStr() + "\t";
}

std::valarray<double> ExponentialCurrentSynapse::GetSynapticState(int pre_neuron)
{
    std::valarray<double> val(1);
    // get average coupling strength
    double Jsum = 0;
    for(unsigned int target=0; target < this->GetNumberOfPostsynapticTargets(pre_neuron); target++){
        Jsum += *(geometry->GetDistributionJ(pre_neuron,target));
    }

    val[0] = Jsum/double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    //val[0] = GetCouplingStrength()*double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    return val;
}

void ExponentialCurrentSynapse::LoadParameters(std::vector<std::string> *input){
    Synapse::LoadParameters(input);

    std::string              name;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);

        if(name.find("exponential_tau") != std::string::npos){
            Tau=(std::stod(values.at(0)));
        }
    }
}


void ExponentialCurrentSynapse::SaveParameters(std::ofstream * stream,std::string id_str){
    Synapse::SaveParameters(stream,id_str);

    *stream << id_str << "exponential_tau\t\t\t\t\t" << std::to_string(Tau) << " seconds\n";
    *stream << "#\t\tThe Post Synaptic Potential decays exponentially through time. The AUC is determined by J and does not depend on tau\n";
}
