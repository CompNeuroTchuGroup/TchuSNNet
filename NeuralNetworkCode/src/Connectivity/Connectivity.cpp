
#include "Connectivity.hpp"
#include "../Synapse/Synapse.hpp"
// THIS IS TO CIRCUMVENT AN INCLUDE LOOP (BAD PRACTICE, BUT THERE IS A NEED FOR CIRCULAR INCLUDES BETWEEN SYNAPSE AND CONNECTIVITY)

Connectivity::Connectivity(Synapse *synapse, GlobalSimInfo *infoGlobal) : infoGlobal{infoGlobal}, synapse{synapse} {
  std::mt19937 generator(synapse->GetSeed());
  SetSeed(generator);
  // SUGGESTION: Use vector of vectors instead of array of vectors
  // https://github.com/saiftyfirst/BP_Demos/blob/master/C%2B%2B/arrayOfVectors_vs_vectorOfVectors.cpp
}

void Connectivity::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
  wParameterStream << idString << "connectivity_type\t\t\t\t\t" << GetTypeStr() << "\n";
}

void Connectivity::SetSeed(std::mt19937 &seedGenerator) {
  // if(!userSeed){
  std::uniform_int_distribution<int> distribution(0, INT32_MAX);
  seed      = distribution(seedGenerator);
  generator = std::mt19937(seed);
  // }
}

void Connectivity::LoadParameters(const std::vector<FileEntry> &connectivityParameters) {
  for (auto &[parameterName, parameterValues] : connectivityParameters) {
    (void)parameterName;
    (void)parameterValues;
    // if(parameterName.find("connectivity_seed") != std::string::npos && infoGlobal->globalSeed != -1){
    //     SetSeed(std::stoi(parameterValues.at(0)));
    // }
  }
}