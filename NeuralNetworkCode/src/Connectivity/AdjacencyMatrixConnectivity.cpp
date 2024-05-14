#include "AdjacencyMatrixConnectivity.hpp"
#include "../Synapse/Synapse.hpp"

AdjacencyMatrixConnectivity::AdjacencyMatrixConnectivity(Synapse *synapse, GlobalSimInfo *infoGlobal) : Connectivity(synapse, infoGlobal) {
  NeuronInt noTargetNeurons = synapse->GetNoTargetNeurons();
  NeuronInt noSourceNeurons = synapse->GetNoSourceNeurons();

  std::string idStr = synapse->GetIdStrWithULine();

  connectivityMatrix.resize(noTargetNeurons);
  for (NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
    connectivityMatrix.at(targetNeuron).resize(noSourceNeurons, 0);
  }

  // for (NeuronInt targetNeuron : std::ranges::views::iota(0,noTargetNeurons)) {
  //     for (NeuronInt sourceNeuron : std::ranges::views::iota(0,noSourceNeurons)) {
  //         connectivityMatrix.at(targetNeuron).at(sourceNeuron) = 0;
  //     }
  // }
  connectionProbFile = infoGlobal->pathToInputFile + "Connectivity_Matrix_" + idStr + ".txt";
}

void AdjacencyMatrixConnectivity::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
  Connectivity::SaveParameters(wParameterStream, idString);
  //*stream << idString << "connectivity_noSourceNeurons " << std::to_string(this->noSourceNeurons) << "\n";
  //*stream << idString << "connectivity_ConnectionProba\t" << std::to_string(this->GetConnectionProbability()) << "\n";
  //*stream << "#" << idString << "connectivity_noSourceNeurons " << std::to_string(this->noSourceNeurons) << "\n";
  wParameterStream << "#\t\t" << IDstringAdjacencyMatrixConnectivity << ": Pre and post population has adjacency matrix.\n";
}

void AdjacencyMatrixConnectivity::LoadParameters(const std::vector<FileEntry> &connectivityParameters) {
  Connectivity::LoadParameters(connectivityParameters);

  for (auto &[parameterName, parameterValues] : connectivityParameters) {
    // if (parameterName.find("seed") != std::string::npos) {
    //     SetSeed(std::stoi(parameterValues.at(0)));
    // }
    if (parameterName.find("noSourceNeurons") != std::string::npos) {
      SetNoSourceNeurons(std::stoi(parameterValues.at(0)));
    }
  }
  GetConnectionWeightsFromFile(connectionProbFile);
}

void AdjacencyMatrixConnectivity::SetNoSourceNeurons(NeuronInt readNoSourceNeurons) {
  if (readNoSourceNeurons < 0) {
    noSourceNeurons = 0;
  } else if (readNoSourceNeurons > synapse->GetNoSourceNeurons()) {
    noSourceNeurons = synapse->GetNoSourceNeurons();
  } else if ((readNoSourceNeurons == synapse->GetNoSourceNeurons()) && (synapse->IsRecurrent())) {
    noSourceNeurons = synapse->GetNoSourceNeurons() - 1;
  } else {
    noSourceNeurons = readNoSourceNeurons;
  }
}

double AdjacencyMatrixConnectivity::GetConnectionProbability() const {
  if (synapse->GetNoSourceNeurons() == 0)
    return 0;
  else
    return static_cast<double>(noSourceNeurons) / static_cast<double>(synapse->GetNoSourceNeurons());
}

void AdjacencyMatrixConnectivity::GetConnectionWeightsFromFile(std::string filepath) {

  std::string readStringLine;

  std::ifstream connectivityStream(filepath);
  std::cout << filepath << std::endl;

  NeuronInt sourceNeuron{};
  NeuronInt noSourceNeurons = synapse->GetNoSourceNeurons();
  NeuronInt noTargetNeurons = synapse->GetNoTargetNeurons();

  std::vector<std::vector<double>> tempMatrix;
  try {
    while (std::getline(connectivityStream, readStringLine)) {
      std::vector<std::string> parameterValues = SplitStringToValues(readStringLine);

      if (sourceNeuron >= noSourceNeurons) {
        throw 4;
      }
      if (static_cast<NeuronInt>(parameterValues.size()) < noTargetNeurons) {
        throw 1;
      }
      if (static_cast<NeuronInt>(parameterValues.size()) > noTargetNeurons) {
        throw 3;
      }

      for (size_t targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
        connectivityMatrix.at(targetNeuron).at(sourceNeuron) = std::stoi(parameterValues.at(targetNeuron));
      }
      sourceNeuron++;
    }
    if (sourceNeuron < noSourceNeurons) {
      throw 2;
    }
  } catch (int exception) {
    connectivityStream.close();
    if (exception == 1) {
      throw "Not enough columns of synaptic connection file!";
    } else if (exception == 2) {
      throw "Not enough rows of synaptic connection file!";
    } else if (exception == 3) {
      throw "Warning: Too many columns of synaptic connection file";
    } else if (exception == 4) {
      throw "Warning: Too many rows of synaptic connection file";
    }
  }
  connectivityStream.close();
}

void AdjacencyMatrixConnectivity::ConnectNeurons() {

  NeuronInt countedSourceNeurons{};
  NeuronInt noTargetNeurons = synapse->GetNoTargetNeurons();
  NeuronInt noSourceNeurons = synapse->GetNoSourceNeurons();

  NeuronInt outputInterval = noTargetNeurons / 10;
  if (outputInterval == 0)
    outputInterval = 1;

  // Iterate through all target neurons
  for (NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
    for (NeuronInt sourceNeuron : std::ranges::views::iota(0, noSourceNeurons)) {

      int connection = connectivityMatrix.at(targetNeuron).at(sourceNeuron);

      if (connection != 0) { // Here we would for loop over the connection index to do multiple synapses
        for (int connectionFormed : std::ranges::views::iota(0, connection)) {
          (void)connectionFormed;
          synapse->AllocateSynapse(targetNeuron, sourceNeuron);
          countedSourceNeurons++;
        }
      }
    }

    // if ((targetNeuron) % outputInterval == 0)
    //     std::cout << 100 * (targetNeuron) / noTargetNeurons << "%\n";
  }
  //     // std::cout << "100%\n";
  SetNoSourceNeurons(countedSourceNeurons / noTargetNeurons);
}
