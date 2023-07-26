#include "PRGSynapseContinuous.hpp"

PRGSynapseContinuous::PRGSynapseContinuous(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal):MongilloSynapseContinuous(targetPop,sourcePop,infoGlobal) {
	l_PRG.resize(GetNoSourceNeurons());
}


void PRGSynapseContinuous::LoadParameters(const std::vector<FileEntry>& synapseParameters){

    MongilloSynapseContinuous::LoadParameters(synapseParameters);

    for(auto& [parameterName, parameterValues] : synapseParameters) {
        if(parameterName.find("prg_M") != std::string::npos){
			M_PRG = std::stod(parameterValues.at(0));
        } else if(parameterName.find("prg_tau_l") != std::string::npos){
            tauL_PRG  = std::stod(parameterValues.at(0));
        } else if(parameterName.find("prg_Delta_tau_f") != std::string::npos){
			deltaTauF_PRG = std::stod(parameterValues.at(0));
        } else if (parameterName.find("prg_Delta_U") != std::string::npos) {
			deltaU_PRG = std::stod(parameterValues.at(0));
		}
    }
}


void PRGSynapseContinuous::SaveParameters(std::ofstream& wParameterStream,std::string idString) const{
    MongilloSynapseContinuous::SaveParameters(wParameterStream,idString);

    wParameterStream << idString << "prg_M\t\t\t\t\t\t" << std::to_string(M_PRG) << " #(probability of 0->1 transition of l per transmitted spike)\n";
    wParameterStream << idString << "prg_tau_l\t\t\t\t\t\t" << std::to_string(tauL_PRG) << " #secs (decay time of l)\n";
    wParameterStream << idString << "prg_Delta_tau_f\t\t\t\t\t" << std::to_string(deltaTauF_PRG) << " #Increase of tauF_Mongillo (as defined in MongilloSynapseContinuous) due to LPA : tauf -> tauF_Mongillo+l*Delta_tau_f\n";
    wParameterStream << idString << "prg_Delta_U\t\t\t\t\t" << std::to_string(deltaU_PRG) << " #Increase of U (as defined in MongilloSynapseContinuous) due to LPA : U -> U+l*Delta_U\n";

    wParameterStream << "#\tsynapticTargets transitions from 1->0 at the rate 1/tau_l (LPA2 unbinding) \n";
    wParameterStream << "#\tsynapticTargets transitions from 0->1 for transmitted spikes (with prob. M) due to ATX-upregulation.\n";
}


std::string PRGSynapseContinuous::GetDataHeader(int dataColumn)
{
	return "#" + std::to_string(dataColumn) + " J_" + GetIdStr() + " (mV)\n"
		+ "#" + std::to_string(dataColumn + 1) + " <x>_" + GetIdStr() + " ( x is given post-spike, as it was afetr last time the presynaptic neuron spiked) \n"
		+ "#" + std::to_string(dataColumn + 2) + " <y>_" + GetIdStr() + " ( y is given post-spike) \n"
		+ "#" + std::to_string(dataColumn + 3) + " <l>_" + GetIdStr() + " ( l is given post-spike) \n"
		+ "#" + std::to_string(dataColumn + 4) + " <pR>_" + GetIdStr() + "\n";
	//+ "#" + std::to_string(data_column+4) + " <xy>_"  + GetIdStr() + "\n"   XY will always return 0 since it is calculated after spiking
}


std::string PRGSynapseContinuous::GetUnhashedDataHeader() const
{
	return "J_" + GetIdStr() + "\t<x>_" + GetIdStr() + "\t<y>_" + GetIdStr() + "\t<l>_" + GetIdStr() + "\t<pR>_" + GetIdStr() + "\t";
}

std::vector<double> PRGSynapseContinuous::GetSynapticState(NeuronInt sourceNeuron) const
{
	std::vector<double> value;
	value.resize(GetNoDataColumns());

    double sumXMinus{ 0 };//the synapses where the spike was transmitted a neurotransmitter count reset to 0
    double sumYMinus{ 0 };//expected value of y before the spike-induced increase in bound calcium probability
    double sumLMinus{ 0 };//expected value of y before the spike-induced increase in bound calcium probability
    double spikeCounter{ 0 }; // XY{0},
    double lMinus{ 0 };
    // double x_minus{ 0 }, y_minus{ 0 }, ; //state of x, y and l before the spike was processed
    NeuronInt noTargetNeurons{this->GetNoTargetedNeurons(sourceNeuron)};
    double Jsum = 0;

    for (NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
        // compute l before spike
        lMinus = (l_PRG.at(sourceNeuron).at(targetNeuron)-M_PRG*y_MongilloC.at(sourceNeuron).at(targetNeuron)*x_MongilloC.at(sourceNeuron).at(targetNeuron))/(1-M_PRG*y_MongilloC.at(sourceNeuron).at(targetNeuron)*x_MongilloC.at(sourceNeuron).at(targetNeuron));
        sumLMinus += lMinus;
        // compute y before spike
        // y_minus = (y_Mongillo.at(sourceNeuron).at(targetNeuron) - u_Mongillo-l_minus*deltaU_PRG)/(1-u_Mongillo-l_minus*deltaU_PRG);
        sumYMinus += (y_MongilloC.at(sourceNeuron).at(targetNeuron) - u_MongilloC-lMinus*deltaU_PRG)/(1-u_MongilloC-lMinus*deltaU_PRG);

        // compute x before spike
        if(y_MongilloC.at(sourceNeuron).at(targetNeuron) == 1){
            // double J_minus = GetCouplingStrength(targetNeuron, sourceNeuron);
            // x_minus = spikeSubmitted.at(sourceNeuron).at(targetNeuron)/(GetCouplingStrength(targetNeuron, sourceNeuron) * y_Mongillo.at(sourceNeuron).at(targetNeuron));
            sumXMinus	+=spikeSubmitted.at(sourceNeuron).at(targetNeuron)/(GetCouplingStrength(targetNeuron,sourceNeuron) * y_MongilloC.at(sourceNeuron).at(targetNeuron));
        }
        else{
            // x_minus = x_Mongillo.at(sourceNeuron).at(targetNeuron)/(1-y_Mongillo.at(sourceNeuron).at(targetNeuron));
            sumXMinus	+=x_MongilloC.at(sourceNeuron).at(targetNeuron)/(1-y_MongilloC.at(sourceNeuron).at(targetNeuron));
        }
        // sumXMinus	+= x_minus;

        // count submitted spikes
        spikeCounter += spikeSubmitted.at(sourceNeuron).at(targetNeuron);
        //XY += x[pre_neuron][i] * y[pre_neuron][i];
    // get average coupling strength
        Jsum += GetDistributionJ(targetNeuron,sourceNeuron);
    }

    value.at(0) = Jsum/static_cast<double>(this->GetNoTargetedNeurons(sourceNeuron));
	//value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));

	value.at(1) = sumXMinus;//the synapses where the spike was transmitted a neurotransmitter count reset to 0
	value.at(2) = sumYMinus;//expected value of y before the spike-induced increase in bound calcium probability
	value.at(3) = sumLMinus;//expected value of y before the spike-induced increase in bound calcium probability
	value.at(4) = spikeCounter;
	//	value[4] = static_cast<double>(XY);
	return value;
}


std::vector<double> PRGSynapseContinuous::AdvectSpikers (NeuronInt spiker) {
	double dtLastSpike = sourcePop->GetTimeSinceLastSpike(spiker); //static_cast<double>(infoGlobal->timeStep - neuronsPre->get_previous_spike_step(spiker))*dtTimestep;

    double exptf;
    double exptd        = exp(-dtLastSpike / tauD_MongilloC);
    double exptl        = exp(-dtLastSpike / tauL_PRG);
    std::vector<double> currents(GetNoTargetedNeurons(spiker), 0.0);
    //double exptf_Dtauf  = exp(-dt_lastSpike / (tauF_Mongillo+Delta_tau_f));
	//double exp_LPA2unbiding = exp(-infoGlobal->dtTimestep / tau_l);
	//double delta_l; //time at which l switches back to 0 (randomly generated at each synapse)

	for (NeuronInt targetNeuron : std::ranges::views::iota(0, static_cast<NeuronInt>(currents.size()))){
        // Compute decay of Calcium based on the value of l after the last spike
        exptf = exp(-dtLastSpike / (tauF_MongilloC + l_PRG.at(spiker).at(targetNeuron)*deltaTauF_PRG));

        // Calcium unbinds with time constant tauF_Mongillo + l* deltatauf in between spikes
        y_MongilloC.at(spiker).at(targetNeuron) = y_MongilloC.at(spiker).at(targetNeuron)*exptf;

        // Neurotransmitter is refilled with time constant tauD in between spikes
        x_MongilloC.at(spiker).at(targetNeuron) = 1 - (1-x_MongilloC.at(spiker).at(targetNeuron))*exptd;

        // LPA unbinds with time constant tau_l in between spikes
        l_PRG.at(spiker).at(targetNeuron) = l_PRG.at(spiker).at(targetNeuron)*exptl;


		//Upon presynaptic spike, the amount of Calcium bound increases
        y_MongilloC.at(spiker).at(targetNeuron) += (u_MongilloC+l_PRG.at(spiker).at(targetNeuron)*deltaU_PRG)*(1-y_MongilloC.at(spiker).at(targetNeuron));

		//Spike transmission
		currents.at(targetNeuron)+=TransmitSpike(targetNeuron, spiker);
	}
    return currents;
}


double PRGSynapseContinuous::TransmitSpike(NeuronInt targetNeuron,NeuronInt spiker){

    // Spike is transmitted and neurotransmitter is released
    double current{MongilloSynapseContinuous::TransmitSpike(targetNeuron, spiker)};

    // LPA2 receptor is filled relative to the strength of the transmitted spike
    l_PRG.at(spiker).at(targetNeuron) += M_PRG*y_MongilloC.at(spiker).at(targetNeuron)*x_MongilloC.at(spiker).at(targetNeuron)*(1-l_PRG.at(spiker).at(targetNeuron));
    return current;
}


void PRGSynapseContinuous::ConnectNeurons() {
	MongilloSynapseContinuous::ConnectNeurons();

	//Assemble a list of source neurons that project onto each postsynaptic neurons
	// this is done by going one by one through the list of postsynaptic neurons for each source neurons
    for (NeuronInt sourceNeuron : std::ranges::views::iota(0, GetNoSourceNeurons())) {
		l_PRG.at(sourceNeuron).resize(targetSpineList.at(sourceNeuron).size(), 0);
        // for (int targetNeuron : range(0,noTargetNeurons)) {
        //     		l_PRG.at(sourceNeuron).at(targetNeuron) = 0;
        // 	}
	}
}
