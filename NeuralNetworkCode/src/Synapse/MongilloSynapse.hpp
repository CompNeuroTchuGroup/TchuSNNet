
#ifndef MONGILLOSYNAPSE_CURRENTBASED
#define MONGILLOSYNAPSE_CURRENTBASED
#include "Synapse.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "../GlobalFunctions.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>

/*struct SpineData_STP {
    bool x = 0;
    bool y = 0;
    bool spike_submitted = 0;
    bool gamma_s = false; //this is only used for PRG Synapses (can be removed for Mongillo only)
};*/

class MongilloSynapse : public Synapse {
protected:

    double u_Mongillo{}; //Calcium binding probability
    double tauF_Mongillo{};  // in units of time steps, calcium unbinding time constant
    double tauD_Mongillo{};  // in units of time steps, neurotransmitter refill time constant

    //std::valarray<std::valarray<SpineData_STP>> synapseData;
    std::vector<std::vector<bool>> x_Mongillo; //Neurotransmitter amount
    std::vector<std::vector<bool>> y_Mongillo;//Calcium binding
    std::vector<std::vector<bool>> spikeSubmitted;
    std::uniform_real_distribution<double> uniformDistribution{0.0, 1.0};


    std::vector<double> AdvectSpikers (NeuronInt spiker) override;
    //void advect_finalize(std::vector<std::vector<double>> * waiting_matrix) override {}
    virtual double TransmitSpike(NeuronInt targetNeuron,NeuronInt spiker);

public:
    MongilloSynapse(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal);
    ~MongilloSynapse() override = default;

    void ConnectNeurons() override;

    //*****************************
    //******* Get Functions *******
    //*****************************
    int                     GetNoDataColumns() const override { return 4; } // J, <y>, <x>, <submitted_spikes>
    std::string             GetDataHeader(int dataColumn) override;
	std::string				GetUnhashedDataHeader() const override;
    std::vector<double>   GetSynapticState(NeuronInt sourceNeuron) const override;
    std::string             GetTypeStr() const override { return IDstringMongilloSynapse; };

    void SaveParameters(std::ofstream& wParameterStream,std::string idString) const override;
    void LoadParameters(const std::vector<FileEntry>& synapseParameters) override;
};


#endif // MONGILLOSYNAPSE_CURRENTBASED
