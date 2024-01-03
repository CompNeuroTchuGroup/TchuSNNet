#ifndef QIF_NEURONPOP_HPP
#define QIF_NEURONPOP_HPP

#include "../NeuronPop.hpp"
#include <iostream>
#include <random>
#include <vector>

// sharpness > 0
// 0 < criticalPotential < thresholdV

class QIFNeuronPop : public NeuronPop {
protected:
  double criticalV{};
  double sharpness{};

public:
  // QIFNeuronPop(double tm, double vr, double vc, double s, double vt, double t);
  QIFNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt neuronID) : NeuronPop(infoGlobal, neuronID) {}

  void        Advect(const std::vector<double> &synaptic_dV);
  void        SaveParameters(std::ofstream &wParameterStream) const;
  std::string GetType() const override { return IDstringQIFNeuron; }
};

#endif // QIFNeuronPop_HPP
