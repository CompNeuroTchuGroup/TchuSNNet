
#ifndef _CORRELATED_POISSON_NEURONPOP_HPP
#define _CORRELATED_POISSON_NEURONPOP_HPP

#include "../../GlobalFunctions.hpp"
#include "../NeuronPop.hpp"
#include <iostream>
#include <random>
#include <string>
#include <vector>

class CorrPoissonNeuronPop : public NeuronPop {
  // Poisson neuron pop where all neurons fire with certain correlation to each other (the correlation is global)
protected:
  double finalTargetRate{}; // target firing rate
  double correlationCoefficient{};

  double totalLambda{}; // probability of firing in one timestep
  double corrLambda{};
  double uncorrLambda{};

  std::vector<NeuronInt>                 neuronIds{}; // Used in the random sampling with fixed firing rate
  std::uniform_real_distribution<double> uniformDistribution;
  std::binomial_distribution<>           uncorrBinomialDistribution;

public:
  CorrPoissonNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt id);
  ~CorrPoissonNeuronPop() override = default;

  void Advect(const std::vector<double> &synaptic_dV);

  std::string GetType() const override { return IDstringCorrelatedPoissonNeuron; }
  void        SaveParameters(std::ofstream &wParameterStream) const;
  void        LoadParameters(const std::vector<FileEntry> &parameters);
  void        PreCalcLambdas();
};

#endif // _POISSON_NEURONPOP_HEADER_
