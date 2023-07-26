#include "EIFNeuronPop.hpp"


void EIFNeuronPop::Advect(const std::vector<double>& synaptic_dV)
{

    //int   totalNeurons = this->GetTotalNeurons();
    //int   neuronPop;
    double dtTimestep = infoGlobal->dtTimestep;

    ClearSpikerVector();
	
    for(NeuronInt neuronID : std::ranges::views::iota(0, noNeurons)){
		//(a) wait for refractory period
        if((previousSpikeDistance.at(neuronID)) <= refractorySteps){
			continue;
		}
		//(b)Advect
		membraneV.at(neuronID) += dtTimestep / membraneVTau * (-(membraneV.at(neuronID) - leakPotential) + dtTimestep / membraneVTau * sharpness * exp((membraneV.at(neuronID) - criticalV) / sharpness)) + synaptic_dV.at(neuronID);
		membraneV.at(neuronID) = fmax(lowerBoundPotential, membraneV.at(neuronID));
		//(c)determine if neuron has spiked
		if (this->membraneV.at(neuronID) > thresholdV) {
			spikerNeurons.push_back(neuronID);
			membraneV.at(neuronID) = resetV;
		}
        //    while(this->potential[i] > thresholdV)
        //        this->potential[i] -= thresholdV-resetPotential;
        //}
    }
	this->AdvectPlasticityModel();
}

void EIFNeuronPop::LoadParameters(const std::vector<FileEntry>& neuronParameters) {

	NeuronPop::LoadParameters(neuronParameters);

	for (auto& [parameterName, parameterValues] : neuronParameters) {
		if (parameterName.find("V_Crit") != std::string::npos) {
			criticalV= std::stod(parameterValues.at(0));
		}
		if (parameterName.find("sharpness") != std::string::npos) {
			sharpness = std::stod(parameterValues.at(0));
		}
		if (parameterName.find("V_lowerbound") != std::string::npos) {
			lowerBoundPotential = std::stod(parameterValues.at(0));
		}
		if (parameterName.find("V_leak") != std::string::npos) {
			leakPotential = std::stod(parameterValues.at(0));
		}
	}
}

void EIFNeuronPop::SaveParameters(std::ofstream& wParameterStream) const{

	std::string idString = "neurons_" + std::to_string(GetId());

	NeuronPop::SaveParameters(wParameterStream);
	wParameterStream << idString + "_V_Crit\t\t\t" << std::to_string(criticalV) << " #mV\n";
	wParameterStream << idString + "_sharpness\t\t\t" << std::to_string(sharpness) << "\n";
	wParameterStream << idString + "_V_lowerbound\t\t\t" << std::to_string(lowerBoundPotential) << " #mV\n";
	wParameterStream << idString + "_V_leak\t\t\t" << std::to_string(leakPotential) << " #mV\n";
    wParameterStream <<  "#\t\tEIF neuron : dV/dt = -(V-Vleak)/tauM + sharpness/tauM * exp((V-V_Crit)/sharpness) + RI/tauM \n";
	wParameterStream <<  "#\t\tVcannot be lower than V_lowerbound";
	wParameterStream <<  "#\t\treset: v = v_reset + (v - v_thresh)\n";
}

/*EIFNeuronPop::EIFNeuronPop(double tm, double vr, double vc, double s, double vt, double t) :
NeuronPop()
{
    this->membraneVDecayTau = tm;
    this->resetPotential = vr;
    this->criticalPotential = vc;
    this->sharpness = s;
    this->thresholdV = vt;
    this->dtTimestep = t;
}*/
