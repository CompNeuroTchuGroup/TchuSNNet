//
// Created by Saif Ahmed on 05.08.21.
//

#ifndef NEURALNETWORK_HETEROPOISSONNEURONPOP_HPP
#define NEURALNETWORK_HETEROPOISSONNEURONPOP_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "HeteroNeuronPop.hpp"
#include "../../GlobalFunctions.hpp"

class HeteroPoissonNeuronPop: public HeteroNeuronPop {

protected:
    double r_target; // target firing rate
    double lambda{};
    bool inputDependant{true};
    int seed;
    std::mt19937 generator;
    std::uniform_real_distribution<double> uniformDistribution;

public:
    HeteroPoissonNeuronPop(GlobalSimInfo * info,int id);
    ~HeteroPoissonNeuronPop() override = default;

    void advect(std::vector<double> * synaptic_dV) override;

    std::string GetType() override;
    void SaveParameters(std::ofstream * stream) override;
    void LoadParameters(std::vector<std::string> *input) override;

};


#endif //NEURALNETWORK_HETEROPOISSONNEURONPOP_HPP
