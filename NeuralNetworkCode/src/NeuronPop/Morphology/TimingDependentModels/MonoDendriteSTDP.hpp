//
// Created by Saif Ahmed on 27.06.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#ifndef _MONODENDRITE_BASE_CLASS_HPP
#define _MONODENDRITE_BASE_CLASS_HPP

#include "../Morphology.hpp"
#include "../SynapseSpines/CoopSynapseSpine.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>

class MonoDendriteSTDP : public Morphology {
protected:
  std::vector<signed long> spikedSpinesId;
  std::vector<bool>        spikedSpines;
  // std::vector<std::pair<signed long, double>> preSpikes;
  // std::vector<double> postSpikes;
  // std::vector<std::pair<signed long, double>> thetaChanges;
  // std::vector<std::pair<signed long, double>> weightChanges;
  // End
  double tauTheta{};   // decay constant of heterosynaptic effects in spines
  double lambdaDist{}; // decay constant of heterosynaptic effects over distance between synapses
  double tauDelay{};   // decay constant of heterosynaptic effects over inter-synapse spike timing difference

  double lastPostSpikeTime{-200};

  double thetaExpDecay{1.0};
  double dendriticLength{}; // this would change in case of more complex dendritic geometry (atm it is a single 1D dendrite)
  double synapticGap{};     //  minimum gap between syanpses along dendrite

  double preFactorLTP{};
  double preFactorLTD{};

  signed long              spineIdGenerator{}; // variable used to allocate new synapses
  double                   nextPos{};
  bool                     stepWeights{};
  std::vector<signed long> weightStepBoundary{};
  std::vector<double>      weightStepValue{};
  signed long              currWightStepId{};

  std::vector<bool> integratePostSpike{};
  std::vector<bool> integratePreSpike{};

  std::vector<CoopSpinePtr> spineDataCoop;

  double initialWeights{1.0};

  bool   decayWeights{false};
  double WeightDecayConstant{1.0};
  double weightExpDecay{};

  WeightNormalization weightNormalization{NOPNormalization};
  double              softMaxMultiplier{2.0};

  double baseLTP{};
  double baseLTD{};

  double incrementLTP{};
  double decrementLTD{};

  double       alpha{};
  double       beta{};
  virtual void NormalizeWeights();
  void         HardNormalize();
  void         SoftMaxNormalize();
  virtual void WeightDecay();
  void         ThetaDecay();
  void         UpdateCooperativity(signed long spineID, signed long neighborId);
  // void pseudoCoop(signed long spineID, signed long neighborId);

  virtual void UpdateLTP(signed long spineID) = 0;
  virtual void UpdateLTD(signed long spineID) = 0;

  virtual double gLTP(double deltaT) const = 0;
  virtual double gLTD(double deltaT) const = 0;

  virtual double aLTP(double theta) const;
  virtual double aLTD(double theta) const;

  double getDistanceEffects(const CoopSynapseSpine *const spineA, const CoopSynapseSpine *const spineB) const;
  double getTimingEffects(const CoopSynapseSpine *const spineA, const CoopSynapseSpine *const spineB) const;

public:
  explicit MonoDendriteSTDP(GlobalSimInfo *infoGlobal);
  ~MonoDendriteSTDP() override = default;

  void Advect() override;
  void Reset() override;
  void RecordPostSpike() override;
  void RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) override;

  void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
  void CheckParameters(const std::vector<FileEntry> &parameters) override;

  virtual BaseSpinePtr AllocateNewSynapse(BranchTargeting &bTargeting) override;

  int    GetMaxGapDelay(int delayPerMicroMeter) override { return dendriticLength * delayPerMicroMeter; };
  double GetSynapticDistanceToSoma(int synapseID) override { return spineDataCoop.at(synapseID)->GetDistToSoma(); };
  // Revirtualization
};

#endif // NEURALNETWORK_MONODENDRITE_H
