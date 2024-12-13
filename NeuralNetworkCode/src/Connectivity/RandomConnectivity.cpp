#include "RandomConnectivity.hpp"
#include "../Synapse/Synapse.hpp"

RandomConnectivity::RandomConnectivity(Synapse *synapse, GlobalSimInfo *infoGlobal) : Connectivity(synapse, infoGlobal) {
}

void RandomConnectivity::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
  Connectivity::SaveParameters(wParameterStream, idString);
  //*stream << idString << "connectivity_noSourceNeurons " << std::to_string(this->noSourceNeurons) << "\n";
  wParameterStream << idString << "connectivity_ConnectionProba\t\t\t" << std::to_string(this->GetConnectionProbability()) << "\n";
  //*stream << "#" << idString << "connectivity_noSourceNeurons " << std::to_string(this->noSourceNeurons) << "\n";
  wParameterStream
      << "#\t\t" << IDstringRandomConnectivity
      << ": Each neuron receives C = connectionProbability*N_p randomly chosen connections from the presynaptic population p (as used by [Brunel (2000)]).\n";
}

void RandomConnectivity::LoadParameters(const std::vector<FileEntry> &connectivityParameters) {
  Connectivity::LoadParameters(connectivityParameters);

  for (auto &[parameterName, parameterValues] : connectivityParameters) {
    // if(parameterName.find("seed") != std::string::npos){
    //     SetSeed(std::stoi(parameterValues.at(0)));
    // }
    if (parameterName.find("ConnectionProba") != std::string::npos || parameterName.find("connectionProbability") != std::string::npos) {
      // std::cout << "number of neurons pre: " << std::to_string(synapse->GetNoNeuronsPre()) << ", ";
      // std::cout << "connection probability value: " << std::to_string(std::stod(parameterValues.at(0))) << ", ";

      // std::cout << "their product: " << std::to_string(synapse->GetNoNeuronsPre()*std::stod(parameterValues.at(0))) << "\n";
      //  SetNoSourceNeurons(static_cast<signed long>((synapse->GetNoSourceNeurons() * std::stod(parameterValues.at(0)))));
      SetNoSourceNeurons(static_cast<signed long>(std::lround(synapse->GetNoSourceNeurons() * std::stod(parameterValues.at(0)))));
    } else if (parameterName.find("noSourceNeurons") != std::string::npos) {
      SetNoSourceNeurons(std::stoi(parameterValues.at(0)));
    }
  }
}

void RandomConnectivity::SetNoSourceNeurons(NeuronInt readNoSourceNeurons) {
  if (readNoSourceNeurons < 0) {
    noSourceNeurons = 0;
  } else if (readNoSourceNeurons > synapse->GetNoSourceNeurons()) {
    noSourceNeurons = synapse->GetNoSourceNeurons();
  } else if ((readNoSourceNeurons == synapse->GetNoSourceNeurons()) && (synapse->IsRecurrent())) {
    // noSourceNeurons = 0;
    noSourceNeurons = synapse->GetNoSourceNeurons() - 1;
  } else {
    noSourceNeurons = readNoSourceNeurons;
  }
}

double RandomConnectivity::GetConnectionProbability() const {
  if (synapse->GetNoSourceNeurons() == 0) {
    return 0;
  } else if ((noSourceNeurons == synapse->GetNoSourceNeurons() - 1) && synapse->IsRecurrent()) {
    // This is necessary to avoid iterative Parameter files from lowering the connection probability in recurrent synapses
    return (static_cast<double>(noSourceNeurons + 1)) / (static_cast<double>(synapse->GetNoSourceNeurons()));
  } else {
    return (static_cast<double>(noSourceNeurons)) / (static_cast<double>(synapse->GetNoSourceNeurons()));
  }
}

void RandomConnectivity::ConnectNeurons() {

  NeuronInt sourceNeuron, countedSourceNeurons;
  NeuronInt noTargetNeurons = synapse->GetNoTargetNeurons();

  NeuronInt outputInterval = noTargetNeurons / 10;
  if (outputInterval == 0) {
    outputInterval = 1;
  }
  std::uniform_int_distribution<NeuronInt> distribution(0, synapse->GetNoSourceNeurons() - 1);

  // every target neuron has a fixed number of sourceNeuron neurons

  // Iterate through all target neurons
  for (NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
    countedSourceNeurons = 0;
    // assign for each target neuron 'noSourceNeurons' sourceNeuron neurons
    while (countedSourceNeurons < noSourceNeurons) {
      sourceNeuron = distribution(generator);
      // Check if sourceNeuron and target neurons are equal
      if ((synapse->IsRecurrent()) && (sourceNeuron == targetNeuron)) {
        continue;
      }
      // Check if target was assigned to the same sourceNeuron already
      if (!synapse->IsSourceVectorEmpty(sourceNeuron) && synapse->IsTargetLastInVector(targetNeuron, sourceNeuron)) {
        continue;
      }
      synapse->AllocateSynapse(targetNeuron, sourceNeuron);
      countedSourceNeurons++;
    }
    // if((targetNeuron)%outputInterval == 0){
    //     std::cout << 100*(targetNeuron)/noTargetNeurons << "%\n";
    // }
  }
  // std::cout << "100%\n";
}
