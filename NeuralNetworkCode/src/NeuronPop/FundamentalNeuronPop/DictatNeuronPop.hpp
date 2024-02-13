//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _DICTAT_NEURON_POP_HPP
#define _DICTAT_NEURON_POP_HPP

#define _CRT_SECURE_NO_WARNINGS

#include "../../GlobalFunctions.hpp"
#include "../NeuronPop.hpp"

// ControlledNeuronPop
// Reader/FileNeuronPop
// DictatNeuronPop WINNER

class DictatNeuronPop : public NeuronPop {
protected:
  std::vector<std::vector<Instruction>> inputInstructions; // A vector of instructions for each neuron of the population
  std::vector<int>                      activeInstructions;
  bool                                  instructionFiringBool{false};
  bool                                  poissonLikeFiringBool{false};

  std::uniform_real_distribution<double> uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);

  std::string   inputFileAddress;
  std::ifstream inputStreamFile; // A stream to read the spikers per time-step
  bool          spikerListFiringBool{false};

public:
  DictatNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt neuronID);
  ~DictatNeuronPop() override = default;

  void LoadParameters(const std::vector<FileEntry> &parameters) override;
  void SaveParameters(std::ofstream &wParameterStream) const override;

  void Advect(const std::vector<double> &synaptic_dV) override;
  // Here I will have to set the spikers depending on the instructions (using modular with fireEveryNsteps)
  // Possibility of instead of instructions, giving the spikers per timestep (probably calculate the raw firing first, then the spikers per timestep)
  // and set them in the Advect function
  void ReadInstructionsFromFile();
  void GenerateRegularSpikersFromInstructions();
  void GeneratePoissonSpikersFromInstructions();
  void ReadSpikersFromFile();

  std::string GetType() const override { return IDstringDictatNeuron; }
  std::string GetInputFileAddress() const { return inputFileAddress; }
  void        CloseInputStream() { inputStreamFile.close(); }
};

#endif