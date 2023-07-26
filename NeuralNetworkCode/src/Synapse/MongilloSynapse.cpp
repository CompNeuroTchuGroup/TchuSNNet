#include "MongilloSynapse.hpp"


MongilloSynapse::MongilloSynapse(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal):Synapse(targetPop,sourcePop,infoGlobal) {
    x_Mongillo.resize(GetNoSourceNeurons());
    y_Mongillo.resize(GetNoSourceNeurons());
    spikeSubmitted.resize(GetNoSourceNeurons());
}


std::vector<double> MongilloSynapse::AdvectSpikers (NeuronInt spiker){
    double lastSpikeTime    = sourcePop->GetTimeSinceLastSpike(spiker); //static_cast<double>(infoGlobal->timeStep - neuronsPre->get_previous_spike_step(spiker))*dtTimestep;
    double expTauF           = exp(-lastSpikeTime/tauF_Mongillo);
    double expTauD           = exp(-lastSpikeTime/tauD_Mongillo);
    NeuronInt noTargets {GetNoTargetedNeurons(spiker)};
    std::vector<double> currents(noTargets);

    std::replace_if(y_Mongillo.at(spiker).begin(), y_Mongillo.at(spiker).end(), [this, expTauF](bool state){
        return state&&(uniformDistribution(generator) < (1.0-expTauF));
    },false);
    std::replace_if(y_Mongillo.at(spiker).begin(), y_Mongillo.at(spiker).end(), [this](bool state){
        return !state&&(uniformDistribution(generator) < u_Mongillo);
    },true);
    std::replace_if(x_Mongillo.at(spiker).begin(), x_Mongillo.at(spiker).end(), [this, expTauD](bool state){
        return !state&&(uniformDistribution(generator) < (1.0-expTauD));
    },true);
    
    //This for loop is PARALLELIZABLE and thread safe, as currents does not change in size and threads do not need to write in the same position ever, but the estimated gains are probably not that good to risk it.
    for(NeuronInt targetNeuron : std::ranges::views::iota(0,static_cast<NeuronInt>(noTargets))) {
        //SpineData_STP * syn = &(synapseData.at(spiker)[target_counter]);
        //In between spikes: unbind Calcium with rate 1/tauF_Mongillo
        // if((y_Mongillo.at(spiker).at(targetNeuron)) && (uniformDistribution(generator) < (1.0-expTauF))){
        //         y_Mongillo.at(spiker).at(targetNeuron) = false;
        // }
        // //In between spikes: refill neurotransmitter with rate 1/tauD
        // if((!x_Mongillo.at(spiker).at(targetNeuron)) && (uniformDistribution(generator) < (1.0-expTauD))){
        //         x_Mongillo.at(spiker).at(targetNeuron) = true;
        // }
        // //Upon presynaptic spike: bind Calcium with probability u
        // if(((!y_Mongillo.at(spiker).at(targetNeuron))) && (uniformDistribution(generator) < u_Mongillo)){
        //         y_Mongillo.at(spiker).at(targetNeuron) = true;
        // }
        //std::cout << "x = " << std::to_string(x.at(spiker)[target_counter]) << "\n";
        //std::cout << "y = " << std::to_string(y.at(spiker)[target_counter]) << "\n";

        //Spike transmission
        if(x_Mongillo.at(spiker).at(targetNeuron) && y_Mongillo.at(spiker).at(targetNeuron)){
            currents.at(targetNeuron)+=TransmitSpike(targetNeuron, spiker);
        } else {
            spikeSubmitted.at(spiker).at(targetNeuron) = false;
        }
    }
    return currents;
}


double MongilloSynapse::TransmitSpike(NeuronInt targetNeuron,NeuronInt spiker){
    // long target                         = geometry->GetTargetList(spikerId)->at(targetId);

    //double J_ij                         = GetCouplingStrength();
    double J_ij                         = GetCouplingStrength(targetNeuron, spiker);
    x_Mongillo.at(spiker).at(targetNeuron)               = false; //Neurotransmitter release
    spikeSubmitted.at(spiker).at(targetNeuron) = true;

    this->cumulatedDV                  += J_ij; //static_cast<double>(spike_submitted.at(spiker).sum())

    return J_ij;
}


void MongilloSynapse::ConnectNeurons() {
    Synapse::ConnectNeurons();
    for (NeuronInt sourceNeuron : std::ranges::views::iota(0,GetNoSourceNeurons())) {
        NeuronInt noTargetNeurons { GetNoTargetedNeurons(sourceNeuron)};
        x_Mongillo.at(sourceNeuron).resize(noTargetNeurons);
        y_Mongillo.at(sourceNeuron).resize(noTargetNeurons);
        spikeSubmitted.at(sourceNeuron).resize(noTargetNeurons);
    }
}


void MongilloSynapse::LoadParameters(const std::vector<FileEntry>& synapseParameters){

    Synapse::LoadParameters(synapseParameters);

    for(auto& [parameterName, parameterValues] : synapseParameters) {
        if(parameterName.find("mongillo_tauF") != std::string::npos){
            tauF_Mongillo  = std::stod(parameterValues.at(0));
        } else if(parameterName.find("mongillo_tauD") != std::string::npos){
            tauD_Mongillo  = std::stod(parameterValues.at(0));
        } else if(parameterName.find("mongillo_U") != std::string::npos){
            u_Mongillo  = std::stod(parameterValues.at(0));
        // } else if(parameterName.find("mongillo_seed") != std::string::npos){
        //     SetSeed(static_cast<int>(std::stod(parameterValues.at(0))));
        }

    }
}

void MongilloSynapse::SaveParameters(std::ofstream& wParameterStream,std::string idString) const{
    Synapse::SaveParameters(wParameterStream,idString);

    wParameterStream << idString << "mongillo_tauF\t\t\t\t\t" << std::to_string(tauF_Mongillo) << " #secs\n";
    wParameterStream << idString << "mongillo_tauD\t\t\t\t\t" << std::to_string(tauD_Mongillo) << " #secs\n";
    wParameterStream << idString << "mongillo_U\t\t\t\t\t\t" << std::to_string(this->u_Mongillo) << "\n";
	// if (infoGlobal->globalSeed == -1) {
	// 	wParameterStream << idString << "mongillo_seed\t\t\t\t\t" << std::to_string(this->seed) << "\n";
	// }
}


std::string MongilloSynapse::GetDataHeader(int dataColumn) {
    return "#" + std::to_string(dataColumn) + " J_" + GetIdStr() + " (mV)\n"
        + "#" + std::to_string(dataColumn+1) + " <x>_"  + GetIdStr() + " ( x is given pre-spike : x = x_postspike + Spikesubmitted ) \n"
        + "#" + std::to_string(dataColumn+2) + " <y>_"  + GetIdStr() + " ( y is given pre-spike : y = (y_postspike - U)/(1-U) ) \n"
        + "#" + std::to_string(dataColumn+3) + " <pR>_"  + GetIdStr() + "\n";
	//+ "#" + std::to_string(data_column+4) + " <xy>_"  + GetIdStr() + "\n"   XY will always return 0 since it is calculated after spiking
}

std::string MongilloSynapse::GetUnhashedDataHeader() const {
	return "J_" + GetIdStr() + "\t<x>_" + GetIdStr() + "\t<y>_" + GetIdStr() + "\t<pR>_" + GetIdStr() + "\t";
}

std::vector<double> MongilloSynapse::GetSynapticState(NeuronInt sourceNeuron) const {
    std::vector<double> value(GetNoDataColumns());

    int xSum{std::accumulate(x_Mongillo.at(sourceNeuron).begin(), x_Mongillo.at(sourceNeuron).end(), 0, [](int accumulator, bool value){return accumulator+value;})};
    int ySum{std::accumulate(y_Mongillo.at(sourceNeuron).begin(), y_Mongillo.at(sourceNeuron).end(), 0, [](int accumulator, bool value){return accumulator+value;})};
    // int XY=0;
    int SpikeSubmitted{std::accumulate(spikeSubmitted.at(sourceNeuron).begin(), spikeSubmitted.at(sourceNeuron).end(), 0, [](int accumulator, bool value){return accumulator+value;})};
    NeuronInt noTargets {GetNoTargetedNeurons(sourceNeuron)};
    double Jsum {};
    for(int targetNeuron : std::ranges::views::iota(0,noTargets)){
        Jsum += (GetDistributionJ(targetNeuron,sourceNeuron));
        //XY += x[pre_neuron][i] * y[pre_neuron][i];
    }

    value.at(0) = Jsum/static_cast<double>(noTargets);
    //value[0]= GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
    value.at(1)= static_cast<double>(xSum+SpikeSubmitted);//the synapses where the spike was transmitted a neurotransmitter count reset to 0
    value.at(2)= static_cast<double>((ySum-u_Mongillo*noTargets)/(1-u_Mongillo));//expected value of y before the spike-induced increase in bound calcium probability
    value.at(3)= static_cast<double>(SpikeSubmitted);
//	value.at(1) = static_cast<double>(XY);
    return value;
}

// void MongilloSynapse::SetSeed(int inputSeed){
//     seed      = inputSeed;
//     generator = std::mt19937(seed);
// }


// void MongilloSynapse::SetSeed(std::shared_ptr<std::mt19937> generator){
//     std::uniform_int_distribution<int> distribution(0,INT32_MAX);
//     SetSeed(distribution(*generator));
//     Synapse::SetSeed(generator);
// }
