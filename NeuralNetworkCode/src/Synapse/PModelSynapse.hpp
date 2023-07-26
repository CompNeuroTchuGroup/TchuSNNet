#ifndef PLASTICITY_MODEL_SYNAPSE_HEADER_
#define PLASTICITY_MODEL_SYNAPSE_HEADER_
#include "Synapse.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "../GlobalFunctions.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>


class PModelSynapse : public Synapse {
protected:
    std::vector<double> AdvectSpikers (NeuronInt spiker) override;
public:
    PModelSynapse(PopPtr targetPop,PopPtr sourcePop,GlobalSimInfo* infoGlobal);

    //*****************************
    //******* Get Functions *******
    //*****************************
    int GetNoDataColumns() const override { return 1;}
    std::string GetDataHeader(int dataColumn) override;
	std::string GetUnhashedDataHeader() const override;
    std::vector<double> GetSynapticState(NeuronInt sourceNeuron) const override;
    std::string GetTypeStr() const override { return IDstringPlasticityModelSynapse;}

    void SaveParameters(std::ofstream& wParameterStream,std::string idString) const override;
    void LoadParameters(const std::vector<FileEntry>& parameters) override;

    void AllocateSynapse(NeuronInt targetNeuron, NeuronInt sourceNeuron) override;
};


#endif // CURRENTSYNAPSE