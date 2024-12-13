//
//  WhiteNoiseStimulus.hpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 02/03/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#ifndef _WHITE_NOISE_STIMULUS_HPP
#define _WHITE_NOISE_STIMULUS_HPP

#include "Stimulus.hpp"
#include <iostream>
#include <random>
#include <vector>
// #include <limits.h>

class WhiteNoiseStimulus : public Stimulus {
protected:
  // signed seed{};
  std::vector<StepStruct> meanCurrent; // indices of all neurons that have emitted a spike in the previous time step
  std::vector<StepStruct> sigmaCurrent;

  StepStruct *currentMeanStep;
  StepStruct *currentSigmaStep;

  std::vector<double> cachedMeans;
  std::vector<double> cachedSigmas;

  std::normal_distribution<double> standardDistribution{0.0, 1.0};
  // std::mt19937 generator;
  void RecalculateMeansCache(StepStruct *meanStep);
  void RecalculateSigmasCache(StepStruct *sigmaStep);
  void PostLoadParameters() override;
  void SetSignalMatrix() override;

public:
  WhiteNoiseStimulus(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &stimulusParameters, GlobalSimInfo *infoGlobal);
  ~WhiteNoiseStimulus() override = default;

  std::string GetType() const override { return IDstringWhitenoiseStimulus; }
  void        Update(std::vector<std::vector<double>> &synaptic_dV) override;

  double GetScaling(PopInt neuronPop) const override;

  void SaveParameters(std::ofstream &wParameterStream) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
};

#endif /* WhiteNoiseStimulus_hpp */
