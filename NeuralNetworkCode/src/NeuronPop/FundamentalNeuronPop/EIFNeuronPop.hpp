#ifndef EIFNeuronPop_HPP
#define EIFNeuronPop_HPP

#include "../NeuronPop.hpp"
#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>


class EIFNeuronPop : public NeuronPop
{
protected:    
    double criticalV{};
    double sharpness{};
    double lowerBoundPotential{};
    double leakPotential{};
public:
    EIFNeuronPop(GlobalSimInfo* infoGlobal,NeuronInt neuronID) : NeuronPop(infoGlobal,neuronID) {}
    
    void Advect(const std::vector<double>& synaptic_dV);
	void LoadParameters(const std::vector<FileEntry>& neuronParameters);
    void SaveParameters(std::ofstream& wParameterStream) const;

    std::string GetType() const override{return IDstringEIFNeuron;}
};

#endif // EIFNeuronPop_HPP