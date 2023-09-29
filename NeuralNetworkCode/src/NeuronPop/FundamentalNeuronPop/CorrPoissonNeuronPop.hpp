
#ifndef _POISSON_NEURONPOP_HEADER_
#define _POISSON_NEURONPOP_HEADER_

#include "../NeuronPop.hpp"
#include "../../GlobalFunctions.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>



class CorrPoissonNeuronPop : public NeuronPop {
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
    CorrPoissonNeuronPop(GlobalSimInfo* infoGlobal,NeuronInt id);
    ~CorrPoissonNeuronPop() override = default;

    void Advect(const std::vector<double>& synaptic_dV);

    std::string GetType() const override{return IDstringPoissonNeuron;}
    void SaveParameters(std::ofstream& wParameterStream) const;
    void LoadParameters(const std::vector<FileEntry>& parameters);
    void PreCalcLambdas();
};


#endif // _POISSON_NEURONPOP_HEADER_
