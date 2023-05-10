#ifndef CURRENTSYNAPSE
#define CURRENTSYNAPSE

#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <valarray>
#include "Synapse.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "../GlobalFunctions.hpp"

class CurrentSynapse : public Synapse
{
protected:
    void advectSpikers (std::vector<double>& currents, long spiker) override;
    //void advect_finalize(std::vector<double> * synaptic_dV, std::vector<std::vector<std::vector<double>>> * waiting_matrix){}
    //void advect_finalize(std::vector<std::vector<double>> * waiting_matrix) override {}

public:
    CurrentSynapse(NeuronPop * postNeurons,NeuronPop * preNeurons,GlobalSimInfo * info);

    //*****************************
    //******* Get Functions *******
    //*****************************
    int GetNumberOfDataColumns() override { return 1;}
    std::string GetDataHeader(int data_column) override;
	std::string GetUnhashedDataHeader() override;
    std::valarray<double> GetSynapticState(int pre_neuron) override;
    const std::string GetTypeStr() override { return str_currentSynapse; }

    void SaveParameters(std::ofstream * stream,std::string id_str) override;
    void LoadParameters(std::vector<std::string> *input) override;
};


#endif // CURRENTSYNAPSE
