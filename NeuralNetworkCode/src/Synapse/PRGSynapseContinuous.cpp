#include "PRGSynapseContinuous.hpp"

PRGSynapseContinuous::PRGSynapseContinuous(NeuronPop * postNeurons,NeuronPop * preNeurons,GlobalSimInfo * info):MongilloSynapseContinuous(postNeurons,preNeurons,info)
{
	M = 0;
	tau_l = 0;
	Delta_tau_f = 0;
	Delta_U = 0;
	l.resize(GetNoNeuronsPre());
}


void PRGSynapseContinuous::LoadParameters(std::vector<std::string> *input){

    MongilloSynapseContinuous::LoadParameters(input);

    std::string              name;
    std::vector<std::string> values;

    for(auto & it : *input) {
        SplitString(&it,&name,&values);

        if(name.find("prg_M") != std::string::npos){
			M = std::stod(values.at(0));
        }
        else if(name.find("prg_tau_l") != std::string::npos){
            tau_l  = std::stod(values.at(0));
        }
        else if(name.find("prg_Delta_tau_f") != std::string::npos){
			Delta_tau_f = std::stod(values.at(0));
        }
		else if (name.find("prg_Delta_U") != std::string::npos) {
			Delta_U = std::stod(values.at(0));
		}
    }
}


void PRGSynapseContinuous::SaveParameters(std::ofstream * stream,std::string id_str){
    MongilloSynapseContinuous::SaveParameters(stream,id_str);

    *stream << id_str << "prg_M\t\t\t\t\t\t\t" << std::to_string(M) << " (probability of 0->1 transition of l per transmitted spike)\n";
    *stream << id_str << "prg_tau_l\t\t\t\t\t\t" << std::to_string(tau_l) << " seconds (decay time of l)\n";
    *stream << id_str << "prg_Delta_tau_f\t\t\t\t" << std::to_string(Delta_tau_f) << " Increase of tau_f (as defined in MongilloSynapseContinuous) due to LPA : tauf -> tau_f+l*Delta_tau_f\n";
    *stream << id_str << "prg_Delta_U\t\t\t\t\t" << std::to_string(Delta_U) << " Increase of U (as defined in MongilloSynapseContinuous) due to LPA : U -> U+l*Delta_U\n";

    *stream << "#\t\tl transitions from 1->0 at the rate 1/tau_l (LPA2 unbinding) \n";
    *stream << "#\t\tl transitions from 0->1 for transmitted spikes (with prob. M) due to ATX-upregulation.\n";
}


std::string PRGSynapseContinuous::GetDataHeader(int data_column)
{
	return "#" + std::to_string(data_column) + " J_" + GetIdStr() + " (mV)\n"
		+ "#" + std::to_string(data_column + 1) + " <x>_" + GetIdStr() + " ( x is given post-spike, as it was afetr last time the presynaptic neuron spiked) \n"
		+ "#" + std::to_string(data_column + 2) + " <y>_" + GetIdStr() + " ( y is given post-spike) \n"
		+ "#" + std::to_string(data_column + 3) + " <l>_" + GetIdStr() + " ( l is given post-spike) \n"
		+ "#" + std::to_string(data_column + 4) + " <pR>_" + GetIdStr() + "\n";
	//+ "#" + std::to_string(data_column+4) + " <xy>_"  + GetIdStr() + "\n"   XY will always return 0 since it is calculated after spiking
}


std::string PRGSynapseContinuous::GetUnhashedDataHeader()
{
	return "J_" + GetIdStr() + "\t<x>_" + GetIdStr() + "\t<y>_" + GetIdStr() + "\t<l>_" + GetIdStr() + "\t<pR>_" + GetIdStr() + "\t";
}

std::valarray<double> PRGSynapseContinuous::GetSynapticState(int pre_neuron)
{
	std::valarray<double> val;
	val.resize(GetNumberOfDataColumns());

    double X{ 0 }, Y{ 0 }, L{ 0 }, SpikeSubmitted{ 0 }; // XY{0},
    double x_minus{ 0 }, y_minus{ 0 }, l_minus{ 0 };
    int N_post{ static_cast<int>(x[pre_neuron].size()) };

    for (int i{ 0 }; i < N_post; i++) {
        // compute l before spike
        l_minus = (l[pre_neuron][i]-M*y[pre_neuron][i]*x[pre_neuron][i])/(1-M*y[pre_neuron][i]*x[pre_neuron][i]);
        L += l_minus;

        // compute y before spike
        y_minus = (y[pre_neuron][i] - u-l_minus*Delta_U)/(1-u-l_minus*Delta_U);
        Y += y_minus;

        // compute x before spike
        if(y[pre_neuron][i] == 1){
            double J_minus = GetCouplingStrength(pre_neuron, i);
            x_minus = spike_submitted[pre_neuron][i]/(J_minus * y[pre_neuron][i]);
        }
        else{
            x_minus = x[pre_neuron][i]/(1-y[pre_neuron][i]);
        }
        X	+= x_minus;

        // count submitted spikes
        SpikeSubmitted += spike_submitted[pre_neuron][i];
        //XY += x[pre_neuron][i] * y[pre_neuron][i];
    }

    // get average coupling strength
    double Jsum = 0;
    for(unsigned int target=0; target < this->GetNumberOfPostsynapticTargets(pre_neuron); target++){
        Jsum += *(geometry->GetDistributionJ(pre_neuron,target));
    }
    val[0] = Jsum/double(this->GetNumberOfPostsynapticTargets(pre_neuron));
	//val[0] = GetCouplingStrength()*double(this->GetNumberOfPostsynapticTargets(pre_neuron));

	val[1] = double(X);//the synapses where the spike was transmitted a neurotransmitter count reset to 0
	val[2] = double(Y);//expected value of y before the spike-induced increase in bound calcium probability
	val[3] = double(L);//expected value of y before the spike-induced increase in bound calcium probability
	val[4] = double(SpikeSubmitted);
	//	val[4] = double(XY);
	return val;
}


void PRGSynapseContinuous::advectSpikers(std::vector<double>& currents, long spiker)
{
	double dt_lastSpike = neuronsPre->GetTimeSinceLastSpike(spiker); //double(info->time_step - neuronsPre->get_previous_spike_step(spiker))*dt;

    double exptf;
    double exptd        = exp(-dt_lastSpike / tau_d);
    double exptl        = exp(-dt_lastSpike / tau_l);

    //double exptf_Dtauf  = exp(-dt_lastSpike / (tau_f+Delta_tau_f));
	//double exp_LPA2unbiding = exp(-info->dt / tau_l);
	//double delta_l; //time at which l switches back to 0 (randomly generated at each synapse)

	for (unsigned long target_counter = 0;target_counter<x[spiker].size();target_counter++)
	{
        // Compute decay of Calcium based on the value of l after the last spike
        exptf = exp(-dt_lastSpike / (tau_f + l[spiker][target_counter]*Delta_tau_f));

        // Calcium unbinds with time constant tau_f + l* deltatauf in between spikes
        y[spiker][target_counter] = y[spiker][target_counter]*exptf;

        // Neurotransmitter is refilled with time constant tauD in between spikes
        x[spiker][target_counter] = 1 - (1-x[spiker][target_counter])*exptd;

        // LPA unbinds with time constant tau_l in between spikes
        l[spiker][target_counter] = l[spiker][target_counter]*exptl;


		//Upon presynaptic spike, the amount of Calcium bound increases
        y[spiker][target_counter] += (u+l[spiker][target_counter]*Delta_U)*(1-y[spiker][target_counter]);

		//Spike transmission
		TransmitSpike(currents, target_counter, spiker);
	}
}


void PRGSynapseContinuous::TransmitSpike(std::vector<double>& currents,long targetId,long spikerId){

    // Spike is transmitted and neurotransmitter is released
    MongilloSynapseContinuous::TransmitSpike(currents, targetId, spikerId);

    // LPA2 receptor is filled relative to the strength of the transmitted spike
    l[spikerId][targetId] += M*y[spikerId][targetId]*x[spikerId][targetId]*(1-l[spikerId][targetId]);
}


void PRGSynapseContinuous::ConnectNeurons()
{
	MongilloSynapseContinuous::ConnectNeurons();

	//Assemble a list of source neurons that project onto each postsynaptic neurons
	// this is done by going one by one through the list of postsynaptic neurons for each source neurons
    for (unsigned long source{ 0 }; source < GetNoNeuronsPre(); source++) {
        long n{ static_cast<long>(geometry->GetTargetList(source)->size()) };
		l[source].resize(n);
        for (int i_n{ 0 }; i_n < n; i_n++) {
            		l[source][i_n] = 0;
        	}
	}
}
