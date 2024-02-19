#ifndef _NO_STIMULUS_HEADER
#define _NO_STIMULUS_HEADER

#include "../GlobalFunctions.hpp"
#include "Stimulus.hpp"

class NoStimulus : Stimulus {
  virtual void SetSignalMatrix();
  virtual void PostLoadParameters();

public:
  NoStimulus(std::shared_ptr<NeuronPopSample> neurons, GlobalSimInfo *infoGlobal);
  virtual ~NoStimulus() = default;

  virtual void        Update(std::vector<std::vector<double>> &synaptic_dV);
  virtual std::string GetType() const { return IDstringNoStimulus; };
  virtual double      GetScaling(PopInt neuronPop) const { return 0.0; };

  virtual void SaveParameters(std::ofstream &wParameterStream) const;
  virtual void LoadParameters(const std::vector<FileEntry> &input);

  // double GetSignalMatrixPoint(PopInt neuronPop, NeuronInt neuron) const { return signalMatrix.at(neuronPop).at(neuron); }
};

#endif
