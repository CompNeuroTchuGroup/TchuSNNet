
#ifndef _POISSON_NEURONPOP_HPP
#define _POISSON_NEURONPOP_HPP

#include "../../GlobalFunctions.hpp"
#include "../NeuronPop.hpp"
#include <iostream>
#include <random>
#include <string>
#include <vector>

class PoissonNeuronPop : public NeuronPop {
protected:
  double targetRate{}; // target firing rate
  double lambda{};     // probability of firing in one timestep
  bool   inputDependant{true};

  std::vector<NeuronInt> neuronIds{}; // Used in the random sampling with fixed firing rate

  std::uniform_real_distribution<double> uniformDistribution;
  std::binomial_distribution<>           binomialDistribution;

public:
  PoissonNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt id);
  ~PoissonNeuronPop() override = default;

  void Advect(const std::vector<double> &synaptic_dV);

  std::string GetType() const override { return IDstringPoissonNeuron; }
  void        SaveParameters(std::ofstream &wParameterStream) const;
  void        LoadParameters(const std::vector<FileEntry> &parameters);
};

#endif // _POISSON_NEURONPOP_HEADER_
