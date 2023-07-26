#ifndef QIFNeuronPop_HPP
#define QIFNeuronPop_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "../NeuronPop.hpp"

// sharpness > 0
// 0 < criticalPotential < thresholdV


class QIFNeuronPop : public NeuronPop
{
protected:
    double criticalV{};
    double sharpness{};
public:
    //QIFNeuronPop(double tm, double vr, double vc, double s, double vt, double t);
    QIFNeuronPop(GlobalSimInfo* infoGlobal,NeuronInt neuronID) : NeuronPop(infoGlobal,neuronID) {}

    void Advect(const std::vector<double>& synaptic_dV);
    void SaveParameters(std::ofstream& wParameterStream) const;
    std::string GetType() const override{return IDstringQIFNeuron;}
};

#endif // QIFNeuronPop_HPP
