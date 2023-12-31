//
//  NeuronPopSample.hpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 20.11.17.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//
#ifndef NeuronPopSample_hpp
#define NeuronPopSample_hpp

#include "GlobalFunctions.hpp"
#include "NeuronPop/NeuronPop.hpp"
#include "NeuronPop/FundamentalNeuronPop/LIFNeuronPop.hpp"
#include "NeuronPop/FundamentalNeuronPop/EIFNeuronPop.hpp"
#include "NeuronPop/FundamentalNeuronPop/QIFNeuronPop.hpp"
#include "NeuronPop/FundamentalNeuronPop/PoissonNeuronPop.hpp"
#include "./NeuronPop/FundamentalNeuronPop/DictatNeuronPop.hpp"

#include <iostream>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>
#include <fstream>
#include <execution>

class NeuronPopSample {
protected:

    GlobalSimInfo* infoGlobal;

    PopInt noPopulations{};
    std::vector<PopPtr> neuronPops;

public:
    NeuronPopSample(std::vector<FileEntry> neuronParameters, GlobalSimInfo* infoGlobal);

    //*******************
    //Get-Functions
    //*******************
    NeuronInt GetTotalNeurons() const { return infoGlobal->totalNeurons; }
    PopInt  GetTotalPopulations() const { return this->noPopulations; }
    NeuronInt GetNeuronsPop(PopInt popId) const { return neuronPops.at(popId)->GetNoNeurons(); }
    PopPtr GetPop(PopInt popId) const { return neuronPops.at(popId); }
    const std::vector<NeuronInt>& GetSpikers(PopInt neuronPop) const { return neuronPops.at(neuronPop)->GetSpikers();}//This function is called in RecordRasterplot()
    const std::vector<NeuronInt>& GetSpikersPrevDt(PopInt neuronPop) const { return neuronPops.at(neuronPop)->GetSpikersPrevDt();}
	double GetXPosition(PopInt neuronPop, NeuronInt neuron) const { return neuronPops.at(neuronPop)->GetXPosition(neuron); }
	double GetYPosition(PopInt neuronPop, NeuronInt neuron) const { return neuronPops.at(neuronPop)->GetYPosition(neuron); }
	double GetPotential(PopInt neuronPop, NeuronInt neuron) const { return neuronPops.at(neuronPop)->GetPotential(neuron); }
    //*******************
    //Set-Functions
    //*******************
    void Advect(std::vector<std::vector<double>>& synaptic_dV);
    void LoadParameters(const std::vector<FileEntry>& neuronParameters);
    void LoadNeuronPop(std::string neuronPopType, PopInt popID, std::vector<FileEntry> neuronParameters);
    void SaveParameters(std::ofstream& wParameterStream) const;
};

#endif /* NeuronPopSample_hpp */
