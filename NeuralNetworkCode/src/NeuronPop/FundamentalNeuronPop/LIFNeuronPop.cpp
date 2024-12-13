
#include "LIFNeuronPop.hpp"

LIFNeuronPop::LIFNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt neuronID) : NeuronPop(infoGlobal, neuronID) {
}

void LIFNeuronPop::Advect(const std::vector<double> &synaptic_dV) {
  ClearSpikerVector();
  // #pragma omp parallel for
  for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
    //(a) wait for refractory period
    if ((previousSpikeDistance.at(neuron)) <= refractorySteps) {
      // std::cout<<(infoGlobal->timeStep-previousSpikeStep.at(neuron));
      continue;
    }

    //(b) Advect
    membraneV.at(neuron) = membraneV.at(neuron) * membraneExpDecay + synaptic_dV.at(neuron);

    //(c) determine if neuron has spiked
    if (membraneV.at(neuron) > thresholdV) {
      spikerNeurons.push_back(neuron);

      // Reset potential
      if (resetType == 0) {
        membraneV.at(neuron) = resetV;
      } else if (resetType == 1) {
        membraneV.at(neuron) -= resetV;
        membraneV.at(neuron) = std::fmod(membraneV.at(neuron), thresholdV - resetV);
        membraneV.at(neuron) += resetV;
        // while(membraneV.at(neuron) > thresholdV){
        //     membraneV.at(neuron) = resetV + (membraneV.at(neuron) - thresholdV);//This could create an infinite loop!
        // }
      }
    }
  }
  this->AdvectPlasticityModel();
}

void LIFNeuronPop::LoadParameters(const std::vector<FileEntry> &neuronParameters) {

  NeuronPop::LoadParameters(neuronParameters);

  for (auto &[parameterName, parameterValues] : neuronParameters) {
    if (parameterName.find("resetType") != std::string::npos) {
      resetType = std::stoi(parameterValues.at(0));
    }
  }
}

void LIFNeuronPop::SaveParameters(std::ofstream &wParameterStream) const {

  std::string idString = "neurons_" + std::to_string(GetId());

  NeuronPop::SaveParameters(wParameterStream);
  wParameterStream << idString + "_resetType\t\t\t\t" << std::to_string(resetType) << "\n";
  wParameterStream << "#\t\tLIF neuron: dV/dt = -V/tauM + RI/tauM \n";
  wParameterStream << "#\t\tresetType 0: v = v_reset\n";
  wParameterStream << "#\t\tresetType 1: v = v_reset + (v - v_thresh) \n";
}
