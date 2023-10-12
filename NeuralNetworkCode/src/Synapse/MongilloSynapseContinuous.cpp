#include "MongilloSynapseContinuous.hpp"

MongilloSynapseContinuous::MongilloSynapseContinuous(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal):Synapse(targetPop,sourcePop,infoGlobal) {
    // SetSeed(0);

    x_MongilloC.resize(GetNoSourceNeurons());
    y_MongilloC.resize(GetNoSourceNeurons());
    spikeSubmitted.resize(GetNoSourceNeurons());

    // uniformDistribution = std::uniform_real_distribution<double>(0.0,1.0);
}


std::vector<double> MongilloSynapseContinuous::AdvectSpikers (NeuronInt spiker){
    double dtLastSpike    = sourcePop->GetTimeSinceLastSpike(spiker); //static_cast<double>(infoGlobal->timeStep - neuronsPre->get_previous_spike_step(spiker))*dtTimestep;
    double exptf           = exp(-dtLastSpike/tauF_MongilloC);
    double exptd           = exp(-dtLastSpike/tauD_MongilloC);
    NeuronInt noTargetNeurons{ this->GetNoTargetedSynapses(spiker) };

    std::vector<double> currents(noTargetNeurons);
    //std::cout << "Synapse " << this->GetIdStr() << "\n";
    std::transform(std::execution::unseq,y_MongilloC.at(spiker).begin(), y_MongilloC.at(spiker).end(), y_MongilloC.at(spiker).begin(),std::bind(std::multiplies<double>(),std::placeholders::_1, exptf));
    std::transform(std::execution::unseq,y_MongilloC.at(spiker).begin(), y_MongilloC.at(spiker).end(), y_MongilloC.at(spiker).begin(),[this](double y_MonC){
        return y_MonC+u_MongilloC*(1-y_MonC);
    });
    std::transform(std::execution::unseq,x_MongilloC.at(spiker).begin(), x_MongilloC.at(spiker).end(), x_MongilloC.at(spiker).begin(),[exptd](double x_MonC){
        return 1 - (1-x_MonC)*exptd;
    });

    for(NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)){
        //std::cout << "spiker: " << std::to_string(spiker) << ", target: " << std::to_string(target_counter) << "\n";

        //In between spikes: unbind Calcium with rate 1/tauF_Mongillo
        //y.at(spiker)[target_counter] = (y.at(spiker)[target_counter]-u_Mongillo)*exptf + u_Mongillo;        // for decay back to u_Mongillo
        // y_Mongillo.at(spiker).at(targetNeuron) *=exptf;   
        // //In between spikes: refill neurotransmitter with rate 1/tauD
		// x_Mongillo.at(spiker).at(targetNeuron) = 1 - (1-x_Mongillo.at(spiker).at(targetNeuron))*exptd;

        // //std::cout << "x_pre = " << std::to_string(x_Mongillo.at(spiker)[target_counter]) << " ";
        // //std::cout << "y_pre = " << std::to_string(y.at(spiker)[target_counter]) << " ";

        // //Upon presynaptic spike: bind Calcium with probability u_Mongillo
        // y_Mongillo.at(spiker).at(targetNeuron) += u_MongilloC*(1-y_Mongillo.at(spiker).at(targetNeuron));

        //Spike transmission
        currents.at(targetNeuron) = TransmitSpike(targetNeuron,spiker);
    }
    return currents;
}


double MongilloSynapseContinuous::TransmitSpike(NeuronInt targetNeuron,NeuronInt spiker){

    // double J_ij                         = GetCouplingStrength(targetNeuron, spiker);
    double current {GetCouplingStrength(targetNeuron, spiker)*x_MongilloC.at(spiker).at(targetNeuron)*y_MongilloC.at(spiker).at(targetNeuron)};

    spikeSubmitted.at(spiker).at(targetNeuron) = current;
    x_MongilloC.at(spiker).at(targetNeuron) *= (1-y_MongilloC.at(spiker).at(targetNeuron));  // neurotransmitter release.

    //std::cout << "x_post = " << std::to_string(x_Mongillo[spikerId][targetId]) << " ";
    //std::cout << "y_post = " << std::to_string(y[spikerId][targetId]) << "\n";
    return current;
}


void MongilloSynapseContinuous::ConnectNeurons() {
    Synapse::ConnectNeurons();

    // test initialization
    std::vector<double> xLocalMax;
    std::vector<double> xLocalMin;
    std::vector<double> yLocalMax;
    std::vector<double> yLocalMin;

    for(NeuronInt sourceNeuron : std::ranges::views::iota(0,GetNoSourceNeurons())){
        NeuronInt totalTargets{ GetNoTargetedSynapses(sourceNeuron) };
        x_MongilloC.at(sourceNeuron).resize(totalTargets, 0.1);
        y_MongilloC.at(sourceNeuron).resize(totalTargets, 0.4);
        spikeSubmitted.at(sourceNeuron).resize(totalTargets);

        // test initialization
        xLocalMax.push_back(*std::max_element(std::begin(x_MongilloC.at(sourceNeuron)),std::end(x_MongilloC.at(sourceNeuron))));
        xLocalMin.push_back(*std::min_element(std::begin(x_MongilloC.at(sourceNeuron)),std::end(x_MongilloC.at(sourceNeuron))));
        yLocalMax.push_back(*std::max_element(std::begin(y_MongilloC.at(sourceNeuron)),std::end(y_MongilloC.at(sourceNeuron))));
        yLocalMin.push_back(*std::min_element(std::begin(y_MongilloC.at(sourceNeuron)),std::end(y_MongilloC.at(sourceNeuron))));
        //std::cout << "max x_Mongillo " << std::to_string(*std::max_element(std::begin(x_Mongillo.at(sourceNeuron)),std::end(x_Mongillo.at(sourceNeuron)))) << "\n";
        //std::cout << "min x_Mongillo " << std::to_string(*std::min_element(std::begin(x_Mongillo.at(sourceNeuron)),std::end(x_Mongillo.at(sourceNeuron)))) << "\n";
        //std::cout << "\n";
    }

    // test initialization
    std::cout << "min x = " << std::to_string(*std::min_element(std::begin(xLocalMin),std::end(xLocalMin))) << " ";
    std::cout << "max x = " << std::to_string(*std::max_element(std::begin(xLocalMax),std::end(xLocalMax))) << " ";
    std::cout << "min y = " << std::to_string(*std::min_element(std::begin(yLocalMin),std::end(yLocalMin))) << " ";
    std::cout << "max y = " << std::to_string(*std::max_element(std::begin(yLocalMax),std::end(yLocalMax))) << " ";
    std::cout << "\n \n";
}


void MongilloSynapseContinuous::LoadParameters(const std::vector<FileEntry>& synapseParameters){

    Synapse::LoadParameters(synapseParameters);

    for(auto& [parameterName, parameterValues] : synapseParameters) {
        if(parameterName.find("mongillo_tauF") != std::string::npos){
            tauF_MongilloC  = std::stod(parameterValues.at(0));
        }
        else if(parameterName.find("mongillo_tauD") != std::string::npos){
            tauD_MongilloC  = std::stod(parameterValues.at(0));
        }
        else if(parameterName.find("mongillo_U") != std::string::npos){
            u_MongilloC  = std::stod(parameterValues.at(0));
        }
        // else if(parameterName.find("mongillo_seed") != std::string::npos){
        //     SetSeed(static_cast<int>(std::stod(parameterValues.at(0))));
        // }
    }
}


void MongilloSynapseContinuous::SaveParameters(std::ofstream& wParameterStream,std::string idString) const{
    Synapse::SaveParameters(wParameterStream,idString);

    wParameterStream << idString << "mongillo_tauF\t\t\t\t\t" << std::to_string(tauF_MongilloC) << " #secs\n";
    wParameterStream << idString << "mongillo_tauD\t\t\t\t\t" << std::to_string(tauD_MongilloC) << " #secs\n";
    wParameterStream << idString << "mongillo_U\t\t\t\t\t\t" << std::to_string(this->u_MongilloC) << "\n";
	// if (infoGlobal->globalSeed == -1) {
	// 	*stream << idString << "mongillo_seed\t\t\t\t\t" << std::to_string(this->seed) << "\n";
	// }
}


std::string MongilloSynapseContinuous::GetDataHeader(int dataColumn) {
    return "#" + std::to_string(dataColumn) + " J_" + GetIdStr() + " (mV)\n"
        + "#" + std::to_string(dataColumn+1) + " <x>_"  + GetIdStr() + " ( x is given pre-spike : x = x_postspike + Spikesubmitted ) \n"
        + "#" + std::to_string(dataColumn+2) + " <y>_"  + GetIdStr() + " ( y is given pre-spike : y = (y_postspike - U)/(1-U) ) \n"
        + "#" + std::to_string(dataColumn+3) + " <pR>_"  + GetIdStr() + "\n";
	//+ "#" + std::to_string(data_column+4) + " <xy>_"  + GetIdStr() + "\n"   XY will always return 0 since it is calculated after spiking
}

std::string MongilloSynapseContinuous::GetUnhashedDataHeader() const {
	return "J_" + GetIdStr() + "\t<x>_" + GetIdStr() + "\t<y>_" + GetIdStr() + "\t<pR>_" + GetIdStr() + "\t";
}

std::vector<double> MongilloSynapseContinuous::GetSynapticState(NeuronInt sourceNeuron) const {
    std::vector<double> outputArray;
    outputArray.resize(GetNoDataColumns());

    double sumXMinus{ 0 };//the synapses where the spike was transmitted a neurotransmitter count reset to 0
    double sumYMinus{ 0 };//expected value of y before the spike-induced increase in bound calcium probability
    // double XY=0;
    double recordedSpikeSubmitted{ 0 };
    //double xMinus{ 0 }, yMinus{0}; //State of x and y before spike was processed
    NeuronInt noTargetNeurons{ this->GetNoTargetedSynapses(sourceNeuron) };
    double Jsum{ 0 };
    
    for(NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)){
        //yMinus = (y_Mongillo.at(sourceNeuron).at(targetNeuron) - u_Mongillo)/(1-u_Mongillo);
        sumYMinus += (y_MongilloC.at(sourceNeuron).at(targetNeuron) - u_MongilloC)/(1-u_MongilloC);
        if(y_MongilloC.at(sourceNeuron).at(targetNeuron) == 1){
            //double J_minus =GetCouplingStrength(sourceNeuron, targetNeuron);
            //xMinus = spikeSubmitted.at(sourceNeuron).at(targetNeuron)/(J_minus * y_Mongillo.at(sourceNeuron).at(targetNeuron));
            sumXMinus	+=spikeSubmitted.at(sourceNeuron).at(targetNeuron)/(GetCouplingStrength(targetNeuron,sourceNeuron) * y_MongilloC.at(sourceNeuron).at(targetNeuron));
        } else{
            // double x_now = x_Mongillo[pre_neuron][i];
            //xMinus = x_Mongillo.at(sourceNeuron).at(targetNeuron)/(1-y_Mongillo.at(sourceNeuron).at(targetNeuron));
            sumXMinus+=x_MongilloC.at(sourceNeuron).at(targetNeuron)/(1-y_MongilloC.at(sourceNeuron).at(targetNeuron));
        }
        // sumXMinus	+= xMinus;
        recordedSpikeSubmitted += spikeSubmitted.at(sourceNeuron).at(targetNeuron);
        //Get average coupling strength
        Jsum += GetDistributionJ(targetNeuron,sourceNeuron);
        //XY += x_Mongillo[pre_neuron][i] * y[pre_neuron][i];
    }

    outputArray.at(0) = Jsum/static_cast<double>(this->GetNoTargetedSynapses(sourceNeuron));
    //value[0]= GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
    outputArray.at(1)= sumXMinus;//the synapses where the spike was transmitted a neurotransmitter count reset to 0
    outputArray.at(2)= sumYMinus;//expected value of y before the spike-induced increase in bound calcium probability
    outputArray.at(3)= recordedSpikeSubmitted;//shouldn't this be an average?
//	value[4] = static_cast<double>(XY);
    return outputArray;
}
