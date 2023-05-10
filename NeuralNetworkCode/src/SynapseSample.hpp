#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <vector>
#include <random>
#include <typeinfo>
#include "GlobalFunctions.hpp"
#include "NeuronPopSample.hpp"
#include "Synapse/CurrentSynapse.hpp"
#include "Synapse/MongilloSynapse.hpp"
#include "Synapse/MongilloSynapseContinuous.hpp"
#include "Synapse/PRGSynapseContinuous.hpp"
#include "Synapse/ExponentialCurrentSynapse.hpp"
#include "Synapse/PowerLawSynapse.hpp"
#include "Synapse/HeteroCurrentSynapse.hpp"

#ifndef SYNAPSESAMPLE_HPP
#define SYNAPSESAMPLE_HPP

//class NeuralNetwork;

class SynapseSample
{
protected:

    int                                  generalSynapseSeed;
    int                                  global_D_max;

    GlobalSimInfo   * info;
    NeuronPopSample * neurons;
    Synapse         *** synapses;

    void LoadParameters(std::vector<std::string> *input);
    void SaveSynapseType(std::string name,std::string type,std::vector<std::string> *input);
public:
    SynapseSample(NeuronPopSample * neurons,std::vector<std::string> *input,GlobalSimInfo *info);
    ~SynapseSample();//Need to make default, doable after moving arrays and pointer resources to std types, smart pointers and references

    void   ConnectNeurons();
    void   WriteConnectivity(std::string filename,int noNeuronsConnectivity);
    void   WriteDistributionD(std::string filename,int noNeuronsDelay);
    void   WriteDistributionJ(std::string filename,int noNeuronsJPot);

    //*******************
    // Get-Functions
    //*******************
    int         GetNumberOfDataColumns();
    int         GetNumberOfDataColumns(int post_population, int pre_population);
    std::string GetDataHeader(int data_column);
	std::string GetUnhashedDataHeader();
    int         GetSeed();
    int         GetMaxD(){return global_D_max;}
    /**
     * Returns the sum of the synaptic state variables of neuron pre_neuron
     * targeting neurons in post_population.
     */
    std::valarray<double>	GetSynapticState(int post_population,int pre_population,int pre_neuron);
	double					GetRecurrentInput(int post_population, int pre_population, int post_neuron);
    double					GetCumulatedDV(int post_population, int pre_population);
    int						GetNumberOfPostsynapticTargets(int post_population,int pre_population,int pre_neuron);

    //**************************************
    // Functions which are redirected to the synapse types
    //**************************************
    //void advect(std::vector<std::vector<double>> * synaptic_dV, std::vector<std::vector<std::vector<double>>> * waiting_matrix);
	void advect(std::vector<std::vector<double>> *  synaptic_dV);
	void Reset();
    void SaveParameters(std::ofstream * stream);

    void Test();
};


#endif // SYNAPSESAMPLE_HPP
