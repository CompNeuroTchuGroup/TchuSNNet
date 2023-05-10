
#ifndef PoissonNeuronPop_HPP
#define PoissonNeuronPop_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "../NeuronPop.hpp"
#include "../../GlobalFunctions.hpp"


class PoissonNeuronPop : public NeuronPop
{
protected:
    double r_target{}; // target firing rate
    double lambda{}; //probability of firing in one timestep
    bool inputDependant{true};
    int seed{};

    std::vector<int> neuronIds{}; //Used in the random sampling with fixed firing rate

    std::mt19937 generator;
    std::uniform_real_distribution<double> uniformDistribution;
    std::binomial_distribution<> binomialDistribution;

public:
    PoissonNeuronPop(GlobalSimInfo * info,int id);
    virtual ~PoissonNeuronPop()=default;

    void advect(std::vector<double> * synaptic_dV);

    std::string GetType(){return str_PoissonNeuron;}
    void SaveParameters(std::ofstream * stream);
    void LoadParameters(std::vector<std::string> *input);
};


#endif // PoissonNeuronPop_HPP
