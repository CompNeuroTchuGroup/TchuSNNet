#include "SynapseSample.hpp"

SynapseSample::SynapseSample(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &synapseParameters, GlobalSimInfo *infoGlobal)
    : infoGlobal{infoGlobal}, neurons{neurons}, totalNeuronPops{neurons->GetTotalPopulations()} {

  synapses      = std::vector<std::vector<SynapsePtr>>(totalNeuronPops);
  synapseStates = std::vector<std::vector<bool>>(totalNeuronPops, std::vector<bool>(totalNeuronPops, false));
  for (std::vector<SynapsePtr> &targetVector : synapses) {
    targetVector.resize(totalNeuronPops);
    for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
      targetVector.at(neuronPop) = nullptr;
    }
  }
  // for (std::vector<synapsePtr>& targetVector : synapses){
  //     for (synapsePtr synapse : targetVector){
  //         //Is this necessary? Not really
  //     }
  // }
  // for(PopInt neuronPop : range(0, totalNeuronPops))
  //     synapses[i] = new Synapse*[totalNeuronPops];

  // for(PopInt neuronPop : range(0, totalNeuronPops)){
  //     for(int j = 0; j < totalNeuronPops; j++){
  //         synapses[i][j] = new CurrentSynapse(neurons->GetPop(i),neurons->GetPop(j),infoGlobal);
  //     }
  // }
  LoadParameters(synapseParameters);
}

void SynapseSample::LoadParameters(const std::vector<FileEntry> &synapseParameters) {

  for (auto &[parameterName, parameterValues] : synapseParameters) {
    if ((parameterName.find("type") != std::string::npos) && ((parameterName.find("connectivity") == std::string::npos))) {
      // I think this line is incorrect. If what we are looking for is the synapse type, this is not it. If we want to see connectivity type (for init
      // reasons) then yes
      SaveSynapseType(parameterName, parameterValues.at(0), synapseParameters);
    }
    // if(parameterName.find("generalSynapseSeed") != std::string::npos){
    //     generalSynapseSeed = std::stoi(parameterValues.at(0));
    // }
  }

  // if(infoGlobal->globalSeed != -1){
  //     std::uniform_int_distribution<int> distribution(0,INT32_MAX);
  //     generalSynapseSeed = distribution(infoGlobal->globalGenerator);
  // }

  // if(generalSynapseSeed >= 0){
  //     //Find totalNeuronPops^2 different seed parameterValues
  //     std::mt19937 generator;
  //     std::uniform_int_distribution<int> distribution(0,INT32_MAX);
  //     generator = std::mt19937(generalSynapseSeed);

  //     for(PopInt neuronPop : range(0, totalNeuronPops)){
  //         for(int j = 0; j < totalNeuronPops; j++){
  //             synapses[i][j]->SetSeed(&generator);
  //         }
  //     }
  // }

  // Check seeds (SEPARATE FUNCTION for each neuronPop)
  for (PopInt popIterator1 : std::ranges::views::iota(0, totalNeuronPops)) {
    for (PopInt popIterator2 : std::ranges::views::iota(0, totalNeuronPops)) {
      for (PopInt popIterator3 : std::ranges::views::iota(0, totalNeuronPops)) {
        for (PopInt popIterator4 : std::ranges::views::iota(0, totalNeuronPops)) {
          if (!GetConnectedState(popIterator1, popIterator2) || !GetConnectedState(popIterator3, popIterator4)) {
            continue;
          }
          if ((synapses.at(popIterator1).at(popIterator2)->GetSeed() == synapses.at(popIterator3).at(popIterator4)->GetSeed()) &&
              (popIterator1 > popIterator3 || popIterator2 > popIterator4)) {
            std::cout << "WARNING! Same Seeds in Connections of Populations " << popIterator1 << "_" << popIterator2 << " and " << popIterator3 << "_"
                      << popIterator4 << std::endl;
          }
        }
      }
    }
  }
  /*
  // Get maximal synaptic delay across all synapse populations
  int curr_max = 0;
  for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
      for(int j = 0; j < totalNeuronPops; j++){
          curr_max = synapses[i][j]->GetMaxD();
          if(global_D_max < curr_max){
              global_D_max = curr_max;
          }
      }
  }
  */
}

void SynapseSample::SaveSynapseType(std::string parameterName, std::string type, const std::vector<FileEntry> &synapseParameters) {

  std::string synapseIdentificator;
  // bool found {false};

  for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
    for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
      synapseIdentificator = std::to_string(targetPop) + "_" + std::to_string(sourcePop);
      if (parameterName.find(synapseIdentificator) != std::string::npos && (parameterName.find("connectivity") == std::string::npos) &&
          (parameterName.find("pmodel") == std::string::npos)) {
        // found = true;
        SetUpSynapse(targetPop, sourcePop, type, FilterStringEntries(synapseParameters, "synapses_" + synapseIdentificator));
        synapseStates.at(targetPop).at(sourcePop) = synapses.at(targetPop).at(sourcePop)->IsConnected();
        return;
      }
    }
  }
}

void SynapseSample::SetUpSynapse(PopInt targetPop, PopInt sourcePop, std::string type, std::vector<FileEntry> synapseParameters) {
  if (type == IDstringCurrentSynapse || type == IDstringHeteroSynapse) {
    synapses.at(targetPop).at(sourcePop) = std::make_unique<CurrentSynapse>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  } else if (type == IDstringMongilloSynapse) {
    synapses.at(targetPop).at(sourcePop) = std::make_unique<MongilloSynapse>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  } else if (type == IDstringMongilloSynapseContinuous) {
    synapses.at(targetPop).at(sourcePop) =
        std::make_unique<MongilloSynapseContinuous>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  } else if (type == IDstringPRGSynapseContinuous) {
    synapses.at(targetPop).at(sourcePop) = std::make_unique<PRGSynapseContinuous>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  } else if (type == IDstringExponentialCurrentSynapse) {
    synapses.at(targetPop).at(sourcePop) =
        std::make_unique<ExponentialCurrentSynapse>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  } else if (type == IDstringPowerLawSynapse) {
    synapses.at(targetPop).at(sourcePop) = std::make_unique<PowerLawSynapse>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  } else if (type == IDstringPlasticityModelSynapse) {
    synapses.at(targetPop).at(sourcePop) = std::make_unique<PModelSynapse>(neurons->GetPop(targetPop), neurons->GetPop(sourcePop), infoGlobal);
  }
  synapses.at(targetPop).at(sourcePop)->LoadParameters(synapseParameters);
}

void SynapseSample::SaveParameters(std::ofstream &wParameterStream) const {

  std::string synapsePrefix;

  wParameterStream << "#*************************************************\n";
  wParameterStream << "#************** Synaptic Parameters **************\n";
  wParameterStream << "#*************************************************\n";

  for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
    for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
      if (this->synapses.at(targetPop).at(sourcePop) != nullptr) {
        synapsePrefix = "synapses_" + std::to_string(targetPop) + "_" + std::to_string(sourcePop) + "_";
        wParameterStream << "#*************************************************\n";
        this->synapses.at(targetPop).at(sourcePop)->SaveParameters(wParameterStream, synapsePrefix);
      } else {
        synapsePrefix = "synapses_" + std::to_string(targetPop) + "_" + std::to_string(sourcePop) + "_";
        wParameterStream << "#*************************************************\n";
        wParameterStream << synapsePrefix << "connected\t\t\t\t\t\t" << std::boolalpha << false << std::noboolalpha << "\n";
      }
    }
  }
}

void SynapseSample::Advect(std::vector<std::vector<double>> &synaptic_dV) {

  // for (Synapse* synapse : connectedSynapses) {
  //             synapse->Advect(synaptic_dV.at(synapse->GetTargetPopID()));
  // }
  std::for_each(PAR, connectedSynapses.begin(), connectedSynapses.end(),
                [&synaptic_dV, this](Synapse *synapse) { synapse->Advect(synaptic_dV.at(synapse->GetTargetPopID()), this->_syndVMutex); });
}

void SynapseSample::Reset() {

  std::for_each(connectedSynapses.begin(), connectedSynapses.end(), [](Synapse *synapse) { synapse->ResetWaitingMatrixEntry(); });
}

void SynapseSample::ConnectNeurons() {
  for (std::vector<SynapsePtr> &targetVector : this->synapses) {
    for (SynapsePtr &synapse : targetVector) {
      if ((synapse != nullptr) &&
          (synapseStates.at(synapse->GetTargetPopID())
               .at(synapse->GetSourcePopID()))) { // from LP the translated condition should be synapse !=nulltpr, and IsConnected()
        connectedSynapses.push_back(synapse.get());
      }
    }
  }
  std::for_each(std::execution::par, connectedSynapses.begin(), connectedSynapses.end(), [](Synapse *synapse) {
    std::cout << "Connecting neuronPop. " << synapse->GetSourcePopID() << " with neuronPop. " << synapse->GetTargetPopID() << "...\n" << std::endl;
    synapse->ConnectNeurons();
    synapse->SetDistributionD();
    synapse->SetDistributionJ();
    synapse->PostConnectNeurons();
    std::cout << "NeuronPop. " << synapse->GetSourcePopID() << " has been connected with neuronPop. " << synapse->GetTargetPopID() << ".\n"
              << std::endl;
  });
}

void SynapseSample::WriteConnectivity(std::string filename, NeuronInt noNeuronsConnectivity) const {

  std::for_each(connectedSynapses.begin(), connectedSynapses.end(), [filename, noNeuronsConnectivity](Synapse *const synapse) {
    synapse->WriteConnectivity(filename + '_' + synapse->GetIdStrWithULine(), noNeuronsConnectivity);
  });
}

void SynapseSample::WriteDistributionD(std::string filename, NeuronInt noNeuronsDelay) const {

  std::for_each(connectedSynapses.begin(), connectedSynapses.end(), [filename, noNeuronsDelay](Synapse *const synapse) {
    synapse->WriteDistributionD(filename + '_' + synapse->GetIdStrWithULine(), noNeuronsDelay);
  });
}

void SynapseSample::WriteDistributionJ(std::string filename, NeuronInt noNeuronsJPot) const {

  std::for_each(connectedSynapses.begin(), connectedSynapses.end(), [filename, noNeuronsJPot](Synapse *const synapse) {
    synapse->WriteDistributionJ(filename + '_' + synapse->GetIdStrWithULine(), noNeuronsJPot);
  });
}

int SynapseSample::GetNoDataColumns(PopInt targetPop, PopInt sourcePop) const {
  return synapses.at(targetPop).at(sourcePop)->GetNoDataColumns();
}

std::string SynapseSample::GetDataHeader(int dataColumn) const {
  std::string dataHeader;
  int         columnCounter{};

  std::for_each(connectedSynapses.begin(), connectedSynapses.end(), [&dataHeader, columnCounter, dataColumn](Synapse *const synapse) mutable {
    dataHeader += synapse->GetDataHeader(dataColumn + columnCounter);
    columnCounter += synapse->GetNoDataColumns();
  });
  return dataHeader;
}

std::string SynapseSample::GetUnhashedDataHeader() const {
  std::string unhashedheader;

  std::for_each(connectedSynapses.begin(), connectedSynapses.end(),
                [&unhashedheader](Synapse *const synapse) { unhashedheader += synapse->GetUnhashedDataHeader(); });
  return unhashedheader;
}

std::vector<double> SynapseSample::GetSynapticState(PopInt targetPop, PopInt sourcePop, NeuronInt sourceNeuron) const {
  return synapses.at(targetPop).at(sourcePop)->GetSynapticState(sourceNeuron);
}

double SynapseSample::GetRecurrentInput(PopInt targetPop, PopInt sourcePop, NeuronInt targetNeuron) const {
  return synapses.at(targetPop).at(sourcePop)->GetRecurrentInput(targetNeuron);
}

double SynapseSample::GetCumulatedDV(PopInt targetPop, PopInt sourcePop) const {
  return synapses.at(targetPop).at(sourcePop)->GetCumulatedDV();
}

NeuronInt SynapseSample::GetNoTargetedNeurons(PopInt targetPop, PopInt sourcePop, NeuronInt sourceNeuron) const {
  return synapses.at(targetPop).at(sourcePop)->GetNoTargetedSynapses(sourceNeuron);
}
