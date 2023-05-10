#ifndef POWERLAWSYNAPSE
#define POWERLAWSYNAPSE

#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <valarray>
#include "Synapse.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "../GlobalFunctions.hpp"


class PowerLawSynapse : public Synapse
{
protected:

	double n;
	int	Naveraging;
	std::valarray<int> spike_count;
	std::valarray<std::valarray<double>> ISI_table;

	void advectSpikers(std::vector<double>& currents, long spiker) override;
	//void advect_finalize(std::vector<std::vector<double>> * waiting_matrix) override {};

public:
	PowerLawSynapse(NeuronPop * postNeurons, NeuronPop * preNeurons, GlobalSimInfo * info);
	~PowerLawSynapse() override = default;


	//*****************************
	//******* Get Functions *******
	//*****************************
	int GetNumberOfDataColumns() override { return 1; }
	std::string GetDataHeader(int data_column) override;
	std::string GetUnhashedDataHeader() override;
	std::valarray<double> GetSynapticState(int pre_neuron) override;
	const std::string GetTypeStr() override { return str_powerlawsynapse; };

	void SaveParameters(std::ofstream * stream, std::string id_str) override;
	void LoadParameters(std::vector<std::string> *input) override;
};


#endif // POWERLAWSYNAPSE
