#include "QIFNeuronPop.hpp"

/*QIFNeuronPop::QIFNeuronPop(double tm, double vr, double vc, double s, double vt, double t) :
NeuronPop()
{
    membraneVDecayTau = tm;
    resetPotential = vr;
    criticalPotential = vc;
    sharpness = s;
    thresholdV = vt;
    dtTimestep = t;
}*/

void QIFNeuronPop::Advect(const std::vector<double>& synaptic_dV){

    //int totalNeurons = GetTotalNeurons();
    // PopInt neuronPop;

    ClearSpikerVector();

    for(NeuronInt neuron : std::ranges::views::iota(0,noNeurons)){
        membraneV.at(neuron) += membraneV.at(neuron) *infoGlobal->dtTimestep/membraneVTau*sharpness*(membraneV.at(neuron)-criticalV)+synaptic_dV.at(neuron);
        if(this->membraneV.at(neuron) > thresholdV){
            spikerNeurons.push_back(neuron);
            while(this->membraneV.at(neuron) > thresholdV)
                this->membraneV.at(neuron)-= thresholdV-resetV;
        }
    }
}


void QIFNeuronPop::SaveParameters(std::ofstream& wParameterStream) const{

    NeuronPop::SaveParameters(wParameterStream);
    wParameterStream <<  "#\t\tQIF neuron (UNDER CONSTRUCTION): dV/dt = V /tauM * sharpness * (V - v_critical) + RI/tauM\n";
    wParameterStream <<  "#\t\treset: v = v_reset + (v - v_thresh)\n";
    wParameterStream <<  "#\t\tUNDER CONSTRUCTION: sharpness/ v_critical not defined\n";

}
