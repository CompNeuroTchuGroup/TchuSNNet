//
//  SpatialPoissonStimulus.h
//  NeuralNetworkCode
//
//  Created by Pierre Ekelmans on 06/04/2018.
//

#ifndef _SPATIAL_POISSON_STIMULUS_HPP
#define _SPATIAL_POISSON_STIMULUS_HPP

#include "Stimulus.hpp"
#include <execution>
#include <iostream>
#include <random>
#include <vector>

class SpatialPoissonStimulus : public Stimulus {
private:
  NeuronInt           noExternalNeurons{1};
  std::vector<double> J_External;
  std::vector<double> lengthScale;
  std::vector<double> peakConnectProb;
  // std::vector<double>         poissonValueTable;    // a table of the poisson distribution for the custom made Poisson generator:

  int isConnectionExact{};

  // int      tableEntries{1000};           // length of poisson_value_table (magic numbers much)
  double inputTau{};
  double attenuation;

  std::vector<double>                              externalPosX;         // position on the x axis
  std::vector<double>                              externalPosY;         // position on the x axis
  std::vector<std::vector<std::vector<NeuronInt>>> externalConnectivity; // Postsynaptic neurons [to external presynaptic neuron i][in population
                                                                         // Pi][list of its postsynaptic neurons in this population]
  std::vector<double> externalCurrents;

  std::vector<StepStruct> stimulusSteps; // This vector contains the following commented ones. First is time step and second is the stimulation
  StepStruct             *currentStep;
  // std::vector<double> stimulusSteps.first;
  // std::vector<double> stimulusSteps.second;

  // std::uniform_int_distribution<int> tablePick;
  std::poisson_distribution<int>           poissonDistr;
  std::uniform_real_distribution<double>   zeroToOneDistr = std::uniform_real_distribution<double>(0, 1);
  std::uniform_int_distribution<NeuronInt> sourceFiringPick;

  void UpdatePoissonTable(); // fills the signal_array
                             // void FillPoissonValueTable(double mu); // fills the poisson_value_table
  void   GenerateConnectivity();
  void   GenerateExactConnectivity();
  void   SetPositions();
  double GetScaling(PopInt neuronPop) const override;
  void   PostLoadParameters() override;
  void   SetSignalMatrix() override;

public:
  SpatialPoissonStimulus(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &stimulusParameters, GlobalSimInfo *infoGlobal);
  ~SpatialPoissonStimulus() override = default;

  //*******************
  // Get-Functions
  //*******************
  double      GetStimulusStepEndTime(int stepNo) const { return static_cast<double>(stimulusSteps.at(stepNo).endTimeStep); }
  double      GetStimulusStep(int stepNo) const { return stimulusSteps.at(stepNo).parameterValues.at(0); };
  int         GetStimulusNoSteps() const { return static_cast<int>(stimulusSteps.size()); }
  std::string GetType() const override { return IDstringSpatialPoissonStimulus; }
  // int     GetTableEntries()          {return tableEntries;}

  //*******************
  // Set-Functions
  //*******************
  void AddStimulusStep(double endTime, double stimStep);

  //*******************************************
  void Update(std::vector<std::vector<double>> &synaptic_dV) override;

  void SaveParameters(std::ofstream &wParameterStream) const override;
  void LoadParameters(const std::vector<FileEntry> &stimulusParameters) override;
  // void    LoadParameters(const std::vector<FileEntry>& parameters,double synapticScaling);
};

#endif /* UncorrelatedPoissonLikeStimulus_h */
