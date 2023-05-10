#include "MongilloSynapseContinuous.hpp"

MongilloSynapseContinuous::MongilloSynapseContinuous(NeuronPop * postNeurons,NeuronPop * preNeurons,GlobalSimInfo * info):Synapse(postNeurons,preNeurons,info)
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


void MongilloSynapseContinuous::advectSpikers (std::vector<double>& currents, long spiker)
{
    double dt_lastSpike    = neuronsPre->GetTimeSinceLastSpike(spiker); //double(info->time_step - neuronsPre->get_previous_spike_step(spiker))*dt;
    double exptf           = exp(-dt_lastSpike/tau_f);
    double exptd           = exp(-dt_lastSpike/tau_d);

    //std::cout << "Synapse " << this->GetIdStr() << "\n";

    for(unsigned long target_counter = 0;target_counter<x[spiker].size();target_counter++)
    {
        //std::cout << "spiker: " << std::to_string(spiker) << ", target: " << std::to_string(target_counter) << "\n";

        //In between spikes: unbind Calcium with rate 1/tau_f
        y[spiker][target_counter] = y[spiker][target_counter]*exptf;                // for decay back to 0
        //y[spiker][target_counter] = (y[spiker][target_counter]-u)*exptf + u;        // for decay back to u

        //In between spikes: refill neurotransmitter with rate 1/tauD
		x[spiker][target_counter] = 1 - (1-x[spiker][target_counter])*exptd;

        //std::cout << "x_pre = " << std::to_string(x[spiker][target_counter]) << " ";
        //std::cout << "y_pre = " << std::to_string(y[spiker][target_counter]) << " ";

        //Upon presynaptic spike: bind Calcium with probability u
        y[spiker][target_counter] = y[spiker][target_counter]+u*(1-y[spiker][target_counter]);

        //Spike transmission
        TransmitSpike(currents, target_counter,spiker);
    }
}


void MongilloSynapseContinuous::TransmitSpike(std::vector<double>& currents, long targetId,long spikerId){

    double J_ij                         = GetCouplingStrength(spikerId, targetId);

    this->cumulatedDV                  += J_ij*x[spikerId][targetId]*y[spikerId][targetId];
    currents[targetId]              += J_ij*x[spikerId][targetId]*y[spikerId][targetId];

    spike_submitted[spikerId][targetId] = J_ij*x[spikerId][targetId]*y[spikerId][targetId];
    x[spikerId][targetId]               = x[spikerId][targetId]*(1-y[spikerId][targetId]);  // neurotransmitter release.

    //std::cout << "x_post = " << std::to_string(x[spikerId][targetId]) << " ";
    //std::cout << "y_post = " << std::to_string(y[spikerId][targetId]) << "\n";
}


void MongilloSynapseContinuous::ConnectNeurons()
{
    Synapse::ConnectNeurons();

    // test initialization
    std::vector<double> x_local_max;
    std::vector<double> x_local_min;
    std::vector<double> y_local_max;
    std::vector<double> y_local_min;

    for(unsigned long source = 0;source < GetNoNeuronsPre();source ++){
        unsigned long n{ static_cast<unsigned long>(geometry->GetTargetList(source)->size()) };
        x[source].resize(n);
        y[source].resize(n);
        spike_submitted[source].resize(n);
        for(unsigned int i_n=0; i_n<n; i_n++){
            x[source][i_n] = 0.1;
            y[source][i_n] = 0.4;
        }

        // test initialization
        x_local_max.push_back(*std::max_element(std::begin(x[source]),std::end(x[source])));
        x_local_min.push_back(*std::min_element(std::begin(x[source]),std::end(x[source])));
        y_local_max.push_back(*std::max_element(std::begin(y[source]),std::end(y[source])));
        y_local_min.push_back(*std::min_element(std::begin(y[source]),std::end(y[source])));
        //std::cout << "max x " << std::to_string(*std::max_element(std::begin(x[source]),std::end(x[source]))) << "\n";
        //std::cout << "min x " << std::to_string(*std::min_element(std::begin(x[source]),std::end(x[source]))) << "\n";
        //std::cout << "\n";
    }

    // test initialization
    std::cout << "min x = " << std::to_string(*std::min_element(std::begin(x_local_min),std::end(x_local_min))) << " ";
    std::cout << "max x = " << std::to_string(*std::max_element(std::begin(x_local_max),std::end(x_local_max))) << " ";
    std::cout << "min y = " << std::to_string(*std::min_element(std::begin(y_local_min),std::end(y_local_min))) << " ";
    std::cout << "max y = " << std::to_string(*std::max_element(std::begin(y_local_max),std::end(y_local_max))) << " ";
    std::cout << "\n \n";
}


void MongilloSynapseContinuous::LoadParameters(std::vector<std::string> *input){

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


void MongilloSynapseContinuous::SaveParameters(std::ofstream * stream,std::string id_str){
    Synapse::SaveParameters(stream,id_str);

    *stream << id_str << "mongillo_tauF\t\t\t\t\t" << std::to_string(tau_f) << " seconds\n";
    *stream << id_str << "mongillo_tauD\t\t\t\t\t" << std::to_string(tau_d) << " seconds\n";
    *stream << id_str << "mongillo_U\t\t\t\t\t\t" << std::to_string(this->u) << "\n";
	if (info->globalSeed == -1) {
		*stream << id_str << "mongillo_seed\t\t\t\t\t" << std::to_string(this->seed) << "\n";
	}
}


std::string MongilloSynapseContinuous::GetDataHeader(int data_column)
{
    return "#" + std::to_string(data_column) + " J_" + GetIdStr() + " (mV)\n"
        + "#" + std::to_string(data_column+1) + " <x>_"  + GetIdStr() + " ( x is given pre-spike : x = x_postspike + Spikesubmitted ) \n"
        + "#" + std::to_string(data_column+2) + " <y>_"  + GetIdStr() + " ( y is given pre-spike : y = (y_postspike - U)/(1-U) ) \n"
        + "#" + std::to_string(data_column+3) + " <pR>_"  + GetIdStr() + "\n";
	//+ "#" + std::to_string(data_column+4) + " <xy>_"  + GetIdStr() + "\n"   XY will always return 0 since it is calculated after spiking
}

std::string MongilloSynapseContinuous::GetUnhashedDataHeader()
{
	return "J_" + GetIdStr() + "\t<x>_" + GetIdStr() + "\t<y>_" + GetIdStr() + "\t<pR>_" + GetIdStr() + "\t";
}

std::valarray<double> MongilloSynapseContinuous::GetSynapticState(int pre_neuron)
{
    std::valarray<double> val;
    val.resize(GetNumberOfDataColumns());

    double X{ 0 };
    double Y{ 0 };
    // double XY=0;
    double SpikeSubmitted{ 0 };

    double x_minus{ 0 }, y_minus{0};
    int N_post{ static_cast<int>(x[pre_neuron].size()) };
    for(int i = 0;i<N_post;i++){
        y_minus = (y[pre_neuron][i] - u)/(1-u);
        Y += y_minus;
        if(y[pre_neuron][i] == 1){
            double J_minus = GetCouplingStrength(pre_neuron, i);
            x_minus = spike_submitted[pre_neuron][i]/(J_minus * y[pre_neuron][i]);
        }
        else{
            // double x_now = x[pre_neuron][i];
            x_minus = x[pre_neuron][i]/(1-y[pre_neuron][i]);
        }
        X	+= x_minus;
        SpikeSubmitted += spike_submitted[pre_neuron][i];
        //XY += x[pre_neuron][i] * y[pre_neuron][i];
    }

    // get average coupling strength
    double Jsum{ 0 };
    for(unsigned int target=0; target < this->GetNumberOfPostsynapticTargets(pre_neuron); target++){
        Jsum += *(geometry->GetDistributionJ(pre_neuron,target));
    }

    val[0] = Jsum/double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    //val[0]= GetCouplingStrength()*double(this->GetNumberOfPostsynapticTargets(pre_neuron));
    val[1]= double(X);//the synapses where the spike was transmitted a neurotransmitter count reset to 0
    val[2]= double(Y);//expected value of y before the spike-induced increase in bound calcium probability
    val[3]= double(SpikeSubmitted);//shouldn't this be an average?
//	val[4] = double(XY);
    return val;
}

void MongilloSynapseContinuous::SetSeed(int s){
    seed      = s;
    generator = std::mt19937(seed);
}


void MongilloSynapseContinuous::SetSeed(std::mt19937 *generator){
    std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    SetSeed(distribution(*generator));
    Synapse::SetSeed(generator);
}
