
#ifndef _CORRELATED_POISSON_NEURONPOP_HEADER_
#define _CORRELATED_POISSON_NEURONPOP_HEADER_

#include "../NeuronPop.hpp"
#include "../../GlobalFunctions.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>



class CorrelatedPoissonNeuronPop : public NeuronPop {
protected:
    double finalTargetRate{}; // target firing rate
    double correlationCoefficient{};

    double totalLambda{}; //probability of firing in one timestep
    double corrLambda{};
    double uncorrLambda{};




    std::vector<NeuronInt> neuronIds{}; //Used in the random sampling with fixed firing rate
    std::uniform_real_distribution<double> uniformDistribution;
    std::binomial_distribution<> uncorrBinomialDistribution;

public:
    CorrelatedPoissonNeuronPop(GlobalSimInfo* infoGlobal,NeuronInt id);
    ~CorrelatedPoissonNeuronPop() override = default;

    void Advect(const std::vector<double>& synaptic_dV);

    std::string GetType() const override{return IDstringCorrelatedPoissonNeuron;}
    void SaveParameters(std::ofstream& wParameterStream) const;
    void LoadParameters(const std::vector<FileEntry>& parameters);
    void PreCalcLambdas();
};


#endif // _POISSON_NEURONPOP_HEADER_
