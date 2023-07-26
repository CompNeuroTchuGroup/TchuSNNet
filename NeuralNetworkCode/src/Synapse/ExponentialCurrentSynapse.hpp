#ifndef EXPONENTIALCURRENTSYNAPSE
#define EXPONENTIALCURRENTSYNAPSE

#include "Synapse.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "../GlobalFunctions.hpp"

#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>


//*****************************************************************************
// Incoming synaptic current to neuron i:
//    ds_i/dtTimestep      = - s_i/tau + sum[delta(t-t_incomingSpikes_i - t_synDelay)]
//    dV_{syn,i}/dtTimestep  = J/tau*s_i(t)
// here, we save outgoing current from each neuron, then sum up respectively
// this can be done due to the linearity of s_i - evolution - equation
//*****************************************************************************

class ExponentialCurrentSynapse : public Synapse {
protected:
    double tauConstant{};
    double expDecayConstant{};
	std::vector<double> AdvectSpikers (NeuronInt spiker) override; // here, only expAddon is advected
	//void advect_finalize(std::vector<std::vector<double>> * waiting_matrix) override {}
    void ResetWaitingMatrixEntry() override;  // synaptic_dV is updated here!
	void ResetcumulatedDV() override;

public:
    ExponentialCurrentSynapse(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal);

    //*****************************
    //******* Get Functions *******
    //*****************************
    int GetNoDataColumns() const override { return 1;}
    std::string GetDataHeader(int dataColumn) override;
	std::string GetUnhashedDataHeader() const override;
    std::vector<double> GetSynapticState(NeuronInt sourceNeuron) const override;
    std::string GetTypeStr() const override { return IDstringExponentialCurrentSynapse; };

    void SaveParameters(std::ofstream& wParameterStream,std::string idString) const override;
    void LoadParameters(const std::vector<FileEntry>& parameters) override;

    //~ExponentialCurrentSynapse(){ delete expAddon;}
};


#endif // EXPONENTIALCURRENTSYNAPSE
