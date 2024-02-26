#include "GlobalStructs.hpp"

DendriticSubRegion::DendriticSubRegion(char regionID, std::vector<int> branchesInRegion):
    regionID { regionID }, branchesInRegion { branchesInRegion } {
  // This is still under work
}


void FileEntry::RemoveCommentsInValues(char commentCharacter) {
  std::vector<std::string> uncommentedValues;
  for (std::string &storedValue : parameterValues) {
    if (storedValue.find(commentCharacter) != std::string::npos) {
      break;
    } else {
      uncommentedValues.push_back(storedValue);
    }
  }
  this->parameterValues = uncommentedValues;
}


Instruction::Instruction(FileEntry inputEntry, double dtTimestep):
    neuronId { std::stoi(inputEntry.parameterValues.at(0)) },
    startTimeStep { std::lround(std::stod(inputEntry.parameterValues.at(1)) / dtTimestep) },
    endTimeStep { std::lround(std::stod(inputEntry.parameterValues.at(2)) / dtTimestep) },
    frequency { std::stod(inputEntry.parameterValues.at(3)) }, firingProbability { frequency * dtTimestep } {
  // Constructor
  if (frequency < std::numeric_limits<double>::epsilon()) {  // Zero comparison to avoid division by zero
    this->off = true;
  } else {
    fireEveryNSteps = std::lround((1 / frequency) / dtTimestep);  // Conversion from frequency to timestep period.
    if (fireEveryNSteps == 0) {                                   // If the frequency is close
      std::cout << "\n"
                << "EXCEPTION: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR NEURON " << std::to_string(neuronId)
                << "\n\n\n"
                << "**********************************";
      throw "EXCEPTION: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR CURRENT DT IN DICTAT INPUT FILE";
    }
  }
}