
#ifndef LIF_NEURONPOP_HPP
#define LIF_NEURONPOP_HPP

#include "../../GlobalFunctions.hpp"
#include "../NeuronPop.hpp"
#include <iostream>
#include <string>
#include <vector>

class LIFNeuronPop : public NeuronPop {
protected:
  int resetType{}; // 0: hard reset, 1: transfer overshoot
public:
  LIFNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt neuronID);
  ~LIFNeuronPop() override = default;

  void        Advect(const std::vector<double> &synaptic_dV) override;
  std::string GetType() const override { return IDstringLIFNeuron; }
  void        SaveParameters(std::ofstream &wParameterStream) const override;
  void        LoadParameters(const std::vector<FileEntry> &parameters) override;
};

#endif // LIFNeuronPop_HPP
