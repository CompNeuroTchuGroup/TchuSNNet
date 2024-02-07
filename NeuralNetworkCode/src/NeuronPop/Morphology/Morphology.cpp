//
// Created by Antoni Bertolin on 14.06.23
//
#include "Morphology.hpp"

Morphology::Morphology(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters)
    : infoGlobal(infoGlobal), totalPostSpikes(0), totalPreSpikes(0) {
  std::uniform_int_distribution<int> distribution(0, INT32_MAX);
  this->seed      = distribution(infoGlobal->globalGenerator);
  this->generator = std::mt19937(this->seed);
  LoadParameters(morphologyParameters);
}

void Morphology::LoadParameters(const std::vector<FileEntry> &morphologyParameters) {

  // checks for correct initialization

  for (auto &[parameterName, parameterValues] : morphologyParameters) {

    if (parameterName.find("seed") != std::string::npos && infoGlobal->globalSeed != -1) {
      this->seed      = std::stoi(parameterValues.at(0));
      this->generator = std::mt19937(this->seed);
      userSeed        = true;
    }
    // include here max and min weights
  }
}

void Morphology::CheckParameters(const std::vector<FileEntry> &parameters) {
  for (auto &[parameterName, parameterValues] : parameters) {

    if (parameterName.find("seed") != std::string::npos && infoGlobal->globalSeed != -1) {
      if (!(seed == std::stoi(parameterValues.at(0)))) {
        throw "Seed was not identical between Synaptic parameters for the plasticity model";
      }
    }
    // include here max and min weights
  }
}

void Morphology::SaveParameters(std::ofstream &wParameterFile, std::string neuronIdentificator) const {
  wParameterFile << "#########From here on is all Heterosynaptic plasticity###########\n";

  wParameterFile << neuronIdentificator << "type\t\t\t" << this->GetType() << "\n";
  if (userSeed) {
    wParameterFile << neuronIdentificator << "seed\t\t\t" << this->seed << "\n";
  }
}

void Morphology::RecordPostSpike() {
  this->totalPostSpikes++;
  this->postSpiked = true;
}

void Morphology::RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) {
  // This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
  // Is there supposed to be a different Inhibitory function?
  this->totalPreSpikes++;
  // STDP Analysis
  // this->preSpikes.emplace_back(spikedSynapseId, this->baseSpineData.at(spikedSynapseId)->lastSpike);
}

std::vector<double> Morphology::GetIndividualSynapticProfile(signed long spineID) const {
  return baseSpineData.at(spineID)->GetIndividualSynapticProfile();
}

std::vector<double> Morphology::GetOverallSynapticProfile() const {
  /*
   * returned array organised as follows:
   * item 1: average synaptic weight
   * item 2: total post spikes
   * item 3: total pre spikes
   * item 4: average plasticity events
   * */
  std::vector<double> dataArray(3);
  size_t              sizeOfSpineData{this->baseSpineData.size()};

  dataArray.at(0) = std::accumulate(this->baseSpineData.begin(), this->baseSpineData.end(), 0.0,
                                    [](double accumulator, const BaseSpinePtr spine) { return accumulator + spine->GetWeightUncoupled(); }) /
                    sizeOfSpineData;
  dataArray.at(1) = this->totalPostSpikes;
  dataArray.at(2) = this->totalPreSpikes;
  return dataArray;
}

std::string Morphology::GetOverallSynapticProfileHeaderInfo() const {
  return std::string("{<average weight>, <total post spikes>, <total pre spikes>}");
}

std::string Morphology::GetIndividualSynapticProfileHeaderInfo() const {
  return baseSpineData.at(0)->GetIndividualSynapticProfileHeaderInfo();
}

double Morphology::GenerateSynapticWeight() {
  double                                 weight{};
  std::uniform_real_distribution<double> distribution(this->minWeight, this->maxWeight);
  if (this->distributeWeights) {
    weight = distribution(generator);
  } else {
    weight = this->initialWeights; // assuming a range of weight between 0 and 2, weight is initialized to midpoint: 1
  }
  // this->weightsSum += weight;
  return weight;
}
/*
void Morphology::triggerStatOut(std::string dirPath) {
    std::ofstream preFile;
    std::ofstream postFile;

    std::ofstream thetasFile;
    std::ofstream weightsFile;

    preFile.open (dirPath + "_pres.dat", std::ofstream::out | std::ofstream::trunc);
    postFile.open (dirPath + "_posts.dat", std::ofstream::out | std::ofstream::trunc);
    thetasFile.open (dirPath + "_thetas.dat", std::ofstream::out | std::ofstream::trunc);
    weightsFile.open (dirPath + "_weights.dat", std::ofstream::out | std::ofstream::trunc);

    for (const auto& line: preSpikes) {
        preFile << line.first << ":" << line.second << std::endl;
    }

    for (const auto& line: postSpikes) {
        postFile << line << std::endl;
    }

    signed long count = 0;
    for (const auto& line: theta_changes) {
        if (count % 2 == 0) {
            thetasFile << "s -- ";
        }
        thetasFile << line.first << ":" << line.second << std::endl;
        count++;
    }

    for (const auto& line: weight_changes) {
        weightsFile << line.first << ":" << line.second << std::endl;
    }

    preFile.close();
    postFile.close();
    thetasFile.close();
    weightsFile.close();

}*/

// void Morphology::printThetasAndWeights() {
//     std::cout << "weights: " << std::endl;
//     for (auto& syn: this->baseSpineData) {
//         std::cout << syn->GetWeight() << ", ";
//     }
//     std::cout << std::endl;

//     std::cout << "thetas: " << std::endl;
//     for (auto& syn: this->baseSpineData) {
//         std::cout << syn->GetTheta() << ", ";
//     }
//     std::cout << std::endl;
// }
