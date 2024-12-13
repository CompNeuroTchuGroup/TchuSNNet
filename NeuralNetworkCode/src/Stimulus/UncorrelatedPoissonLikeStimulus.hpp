//
//  UncorrelatedPoissonLikeStimulus.h
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 23/01/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#ifndef _UNCORRELATED_POISSONLIKE_STIMULUS_HPP
#define _UNCORRELATED_POISSONLIKE_STIMULUS_HPP

#include "Stimulus.hpp"
#include <iostream>
#include <random>
#include <vector>

/* class UncorrelatedPoissonLikeStimulus is a public Stimulus
 * It operates like a virtual external neuron distribution
 * that fires poisson like and that projects to all totalNeurons neurons
 * in the network. Every neuron in the network receives input from
 * connectedNeurons uncorrelated external neurons.
 * The firing rate of this external population is piecewise constant during
 * the specified time intervals.
 * The number of input spikes for all network neurons for the next time step
 * is evaluated at once and stored in the double * signal_array.
 * Note that not a poisson generator from the standard library is used but
 * a custom made that is less precise due to discretization but much faster for
 * higher firing rates.
 */

class UncorrelatedPoissonLikeStimulus : public Stimulus {
private:
  NeuronInt           noExternalNeurons{1};
  std::vector<double> J_External;
  std::vector<double> externalCurrents;
  // std::vector<double> poissonValueTable;    // a table of the poisson distribution for the custom made poisson generator:
  // int      tableEntries;           // length of poissonValueTable
  // int      seed;

  std::vector<StepStruct> stimulusSteps; // This vector contains the following commented ones. First is time step and second is the stimulation
  StepStruct             *currentStep;
  // std::vector<double> nextStimTimeStep;
  // std::vector<double>        nextStimStep;

  std::poisson_distribution<int> poissonDistr;
  // std::mt19937 generator;
  // std::uniform_int_distribution<int> distribution;

  void UpdatePoissonTable(); // fills the signal_array
  // inline void FillPoissonValueTable(double mu); // fills the poissonValueTable

  double GetScaling(PopInt neuronPop) const override;
  void   PostLoadParameters() override;
  void   SetSignalMatrix() override;

public:
  UncorrelatedPoissonLikeStimulus(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &stimulusParameters, GlobalSimInfo *infoGlobal);
  ~UncorrelatedPoissonLikeStimulus() override = default;

  //*******************
  // Get-Functions
  //*******************
  double      GetStimulusStepEndTime(int stepNo) const { return stimulusSteps.at(stepNo).endTimeStep; }
  double      GetStimulusStep(int stepNo) const { return stimulusSteps.at(stepNo).parameterValues.at(0); };
  long        GetStimulusNoSteps() const { return static_cast<long>(stimulusSteps.size()); }
  std::string GetType() const override { return IDstringUncorrelatedStimulus; }
  // int     GetTableEntries()          {return tableEntries;}

  //*******************
  // Set-Functions
  //*******************
  // void SetSeed(int seed){this->seed = seed; generator = std::mt19937(seed);};
  void AddStimulusStep(double endTime, double stimStep);

  //*******************************************
  void Update(std::vector<std::vector<double>> &synaptic_dV) override;

  void SaveParameters(std::ofstream &wParameterStream) const override;
  void LoadParameters(const std::vector<FileEntry> &stimulusParameters) override;
  // void    LoadParameters(const std::vector<FileEntry>& parameters,double synapticScaling);
};

#endif /* UncorrelatedPoissonLikeStimulus_h */
