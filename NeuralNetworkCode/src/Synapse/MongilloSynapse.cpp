#include "MongilloSynapse.hpp"


MongilloSynapse::MongilloSynapse(NeuronPop * postNeurons,NeuronPop * preNeurons,GlobalSimInfo * info):Synapse(postNeurons,preNeurons,info)
{
    u                 = 0;
    tau_f             = 0;
    tau_d             = 0;
    SetSeed(0);

    x.resize(GetNoNeuronsPre());
    y.resize(GetNoNeuronsPre());
    spike_submitted.resize(GetNoNeuronsPre());

    uniformDistribution = std::uniform_real_distribution<double>(0.0,1.0);
}


void MongilloSynapse::advectSpikers (std::vector<double>& currents, long spiker)
{
    double dt_lastSpike    = neuronsPre->GetTimeSinceLastSpike(spiker); //double(info->time_step - neuronsPre->get_previous_spike_step(spiker))*dt;
    double exptf           = exp(-dt_lastSpike/tau_f);
    double exptd           = exp(-dt_lastSpike/tau_d);

    for(std::size_t target_counter = 0;target_counter<x[spiker].size();target_counter++)
    {
        //SynapseData_STP * syn = &(synapseData[spiker][target_counter]);
        //In between spikes: unbind Calcium with rate 1/tau_f
        if((y[spiker][target_counter]) && (uniformDistribution(generator) < (1.0-exptf)))
                y[spiker][target_counter] = false;

        //In between spikes: refill neurotransmitter with rate 1/tauD
        if((!x[spiker][target_counter]) && (uniformDistribution(generator) < (1.0-exptd)))
                x[spiker][target_counter] = true;

        //Upon presynaptic spike: bind Calcium with probability u
        if(((!y[spiker][target_counter])) && (uniformDistribution(generator) < u))
                y[spiker][target_counter] = true;

        //std::cout << "x = " << std::to_string(x[spiker][target_counter]) << "\n";
        //std::cout << "y = " << std::to_string(y[spiker][target_counter]) << "\n";

        //Spike transmission
        if(x[spiker][target_counter] && y[spiker][target_counter])
            TransmitSpike(currents, static_cast<long>(target_counter), spiker);
        else
            spike_submitted[spiker][target_counter] = false;
    }
}


void MongilloSynapse::TransmitSpike(std::vector<double>& currents, long targetId,long spikerId){
    // long target                         = geometry->GetTargetList(spikerId)->at(targetId);

    //double J_ij                         = GetCouplingStrength();
    double J_ij                         = GetCouplingStrength(spikerId, targetId);
    x[spikerId][targetId]               = false; //Neurotransmitter release
    spike_submitted[spikerId][targetId] = true;

    this->cumulatedDV                  += J_ij; //double(spike_submitted[spiker].sum())

    currents[targetId] += J_ij;
}


void MongilloSynapse::ConnectNeurons()
{
    Synapse::ConnectNeurons();
    for (unsigned long source = 0; source < GetNoNeuronsPre(); source++) {
        long n { static_cast<long>(geometry->GetTargetList(source)->size())};
        x[source].resize(n);
        y[source].resize(n);
        spike_submitted[source].resize(n);
    }
}


void MongilloSynapse::LoadParameters(std::vector<std::string> *input){

    Synapse::LoadParameters(input);

    std::string              name;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);

        if(name.find("mongillo_tauF") != std::string::npos){
            tau_f  = std::stod(values.at(0));
        }
        else if(name.find("mongillo_tauD") != std::string::npos){
            tau_d  = std::stod(values.at(0));
        }
        else if(name.find("mongillo_U") != std::string::npos){
            u  = std::stod(values.at(0));
        }
        else if(name.find("mongillo_seed") != std::string::npos){
            SetSeed(static_cast<int>(std::stod(values.at(0))));
        }

    }
}

void MongilloSynapse::SaveParameters(std::ofstream * stream,std::string id_str){
    Synapse::SaveParameters(stream,id_str);

    *stream << id_str << "mongillo_tauF\t\t\t\t\t" << std::to_string(tau_f) << " seconds\n";
    *stream << id_str << "mongillo_tauD\t\t\t\t\t" << std::to_string(tau_d) << " seconds\n";
    *stream << id_str << "mongillo_U\t\t\t\t\t\t" << std::to_string(this->u) << "\n";
	if (info->globalSeed == -1) {
		*stream << id_str << "mongillo_seed\t\t\t\t\t" << std::to_string(this->seed) << "\n";
	}
}


std::string MongilloSynapse::GetDataHeader(int data_column)
{
    return "#" + std::to_string(data_column) + " J_" + GetIdStr() + " (mV)\n"
        + "#" + std::to_string(data_column+1) + " <x>_"  + GetIdStr() + " ( x is given pre-spike : x = x_postspike + Spikesubmitted ) \n"
        + "#" + std::to_string(data_column+2) + " <y>_"  + GetIdStr() + " ( y is given pre-spike : y = (y_postspike - U)/(1-U) ) \n"
        + "#" + std::to_string(data_column+3) + " <pR>_"  + GetIdStr() + "\n";
	//+ "#" + std::to_string(data_column+4) + " <xy>_"  + GetIdStr() + "\n"   XY will always return 0 since it is calculated after spiking
}

std::string MongilloSynapse::GetUnhashedDataHeader()
{
	return "J_" + GetIdStr() + "\t<x>_" + GetIdStr() + "\t<y>_" + GetIdStr() + "\t<pR>_" + GetIdStr() + "\t";
}

std::valarray<double> MongilloSynapse::GetSynapticState(int pre_neuron)
{
    std::valarray<double> val;
    val.resize(GetNumberOfDataColumns());

    int X=0;
    int Y=0;
    // int XY=0;
    int SpikeSubmitted=0;
    int N_post{ static_cast<int>(x[pre_neuron].size()) };

    for(int i = 0;i<N_post;i++){
        X	+= x[pre_neuron][i];
        Y	+= y[pre_neuron][i];
        SpikeSubmitted += spike_submitted[pre_neuron][i];
        //XY += x[pre_neuron][i] * y[pre_neuron][i];
    }

    // get average coupling strength
    double Jsum = 0;
    for(unsigned int target=0; target < this->GetNumberOfPostsynapticTargets(pre_neuron); target++){
        Jsum += *(geometry->GetDistributionJ(pre_neuron,target));
    }

    val[0] = Jsum/double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    //val[0]= GetCouplingStrength()*double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    val[1]= double(X+SpikeSubmitted);//the synapses where the spike was transmitted a neurotransmitter count reset to 0
    val[2]= double((Y-u*N_post)/(1-u));//expected value of y before the spike-induced increase in bound calcium probability
    val[3]= double(SpikeSubmitted);
//	val[4] = double(XY);
    return val;
}

void MongilloSynapse::SetSeed(int s){
    seed      = s;
    generator = std::mt19937(seed);
}


void MongilloSynapse::SetSeed(std::mt19937 *generator){
    std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    SetSeed(distribution(*generator));
    Synapse::SetSeed(generator);
}
