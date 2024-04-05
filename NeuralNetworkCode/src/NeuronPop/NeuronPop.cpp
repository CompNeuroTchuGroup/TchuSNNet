
#include "NeuronPop.hpp"

NeuronPop::NeuronPop(GlobalSimInfo *globalInfo, PopInt inputID): infoGlobal { globalInfo }, identifier { inputID } {
  std::uniform_int_distribution<int> distribution(0, INT_MAX);
  seed = distribution(infoGlobal->globalGenerator);
  // seedInitialPreviousSpike  = distribution(infoGlobal->globalGenerator);
  // seedInitialPotentials = distribution(infoGlobal->globalGenerator);
  generator = std::mt19937(seed);
}

SynInt NeuronPop::GetNoSynapses() const {
  return std::reduce(morphology.begin(), morphology.end(), static_cast<SynInt>(0),
                     [](SynInt accumulator, const std::unique_ptr<Morphology> &morpho) { return accumulator + morpho->GetNoSynapses(); });
}

void NeuronPop::SetNeurons() {
  // std::mt19937 generatorPreviousSpikes(seedInitialPreviousSpike);

  previousSpikeDistance.resize(noNeurons);
  membraneV.resize(noNeurons);

  // std::mt19937 generatorPotentials(seedInitialPotentials);
  // std::uniform_real_distribution<double> uniformDistribution (0.0,1.0); //thresholdV instead of 1.0
  std::uniform_real_distribution<double> potentialDistribution(resetV, thresholdV);
  std::uniform_real_distribution<double> timeDistribution(
    infoGlobal->dtTimestep, (1 / prevMinFrequency));  // The expectation of exponential is the inverse of the lambda. By giving avg freq as
                                                      // lambda, we get avg ISI as expectation (in seconds).

  for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
    previousSpikeDistance.at(neuron) = static_cast<TStepInt>(timeDistribution(generator) / infoGlobal->dtTimestep);
    membraneV.at(neuron)             = potentialDistribution(generator);
  }
}

void NeuronPop::SetPosition() {
  if (infoGlobal->dimensions == 0) {
    return;
  }
  // std::mt19937 generator(seedInitialPotentials);
  std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
  xAxisPosition.resize(noNeurons);
  yAxisPosition.resize(noNeurons);
  double yAxisLength = infoGlobal->yAxisLength;
  double xAxisLength = infoGlobal->xAxisLength;
  int    rows;
  int    columns;
  int    mod;

  // for (int i = 0; i < nos; i++) {
  //	xAxisPosition[i] = uniformDistribution(generator) * infoGlobal->xAxisLength;
  //	yAxisPosition[i] = uniformDistribution(generator) * infoGlobal->yAxisLength;
  // }

  if (infoGlobal->dimensions == 2) {
    rows    = static_cast<int>(ceil(sqrt(noNeurons)));
    columns = noNeurons / rows;
    mod     = noNeurons - rows * columns;  // the difference is added by including one extra column in each of the first mod rows

    for (NeuronInt neuron : std::ranges::views::iota(0, mod * (columns + 1))) {
      // This loop only runs if it is not a perfect square, shifts the available grid into one slightly tighter
      xAxisPosition.at(neuron) = (xAxisLength / (columns + 1)) * (neuron % (columns + 1));  // Assign the position of the column (x axis)
      yAxisPosition.at(neuron) = (yAxisLength / rows) * floor(neuron / (columns + 1));      // The same thing for y
    }
    for (NeuronInt neuron : std::ranges::views::iota(mod * (columns + 1), noNeurons)) {
      // This loop runs from the endpoint of the last loop to the total number of neurons
      xAxisPosition.at(neuron) = (xAxisLength / (columns)) * ((neuron - mod) % columns);
      yAxisPosition.at(neuron) = (yAxisLength / rows) * floor((neuron - mod) / (columns));
    }
  } else if (infoGlobal->dimensions == 1) {
    for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
      xAxisPosition.at(neuron) = neuron * xAxisLength / noNeurons;
    }
  }
}

void NeuronPop::ClearSpikerVector() {
  std::for_each(spikerNeurons.begin(), spikerNeurons.end(),
                [this](NeuronInt previousSpiker) { previousSpikeDistance.at(previousSpiker) = 0; });
  // v1
  //  spikerNeuronsPrevdt=std::move(spikerNeurons);
  //  spikerNeurons.clear();
  // v2. This is done to preserve capacity of spiker vector. It would be faster to just have pointers and swap those, but that would add
  // dereferences everywhere (tested and faster in MVSC)
  std::swap(spikerNeurons, spikerNeuronsPrevdt);
  spikerNeurons.clear();
  std::transform(PAR_UNSEQ, previousSpikeDistance.begin(), previousSpikeDistance.end(), previousSpikeDistance.begin(),
                 std::bind(std::plus<TStepInt>(), std::placeholders::_1, 1));
}

void NeuronPop::AdvectPlasticityModel() {
  if (this->hasPlasticity) {
    std::for_each(PAR_UNSEQ, spikerNeurons.begin(), spikerNeurons.end(), [this](NeuronInt spiker) { RecordPostSpike(spiker); });
    std::for_each(PAR_UNSEQ, morphology.begin(), morphology.end(),
                  [](std::unique_ptr<Morphology> &singleMorphology) { singleMorphology->Advect(); });
  }
}

// void NeuronPop::LoadParameters(const std::vector<FileEntry>& neuronParameters, signed long noNeurons) {
//     SetNeurons(noNeurons);
//     LoadParameters(neuronParameters);
// }

void NeuronPop::LoadParameters(const std::vector<FileEntry> &neuronParameters) {
  for (auto &[parameterName, parameterValues] : neuronParameters) {
    // This is what is called a "structured binding", it assigns the first value to the first variable, etc for the whole container.
    /*if(parameterName.find("Ni") != std::string::npos){
        totalPopulations = (int)parameterValues.size();
        //spiker              = new std::vector<int>[totalPopulations];
        //neuronsInPopulation = new int[totalPopulations];
        //for(int i = 0;i<totalPopulations;i++)
        //    neuronsInPopulation[i] = std::stoi(parameterValues.at(i));
    }*/
    if (parameterName.find("tauM") != std::string::npos) {
      this->membraneVTau     = (std::stod(parameterValues.at(0)));
      this->membraneExpDecay = exp(-infoGlobal->dtTimestep / membraneVTau);
    } else if (parameterName.find("vReset") != std::string::npos) {
      this->resetV = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("vThresh") != std::string::npos) {
      this->thresholdV = (std::stod(parameterValues.at(0)));
    } else if (parameterName.find("refractoryTime") != std::string::npos) {
      this->refractorySteps = std::lround((std::stod(parameterValues.at(0))) / infoGlobal->dtTimestep);
      //} else if(parameterName.find("r_target") != std::string::npos){
      //   this->targetRate = (std::stoi(parameterValues.at(0)));
    } else if (parameterName.find("seed") != std::string::npos) {
      this->userSeeds = true;
      this->seed      = (std::stoi(parameterValues.at(0)));
    } else if (parameterName.find("prevMinFrequency") != std::string::npos) {
      this->prevMinFrequency    = (std::stod(parameterValues.at(0)));
      this->definedPrevMeanFreq = true;
      // } else if(parameterName.find("seedInitialPotentials") != std::string::npos){
      //     userSeeds=true;
      //     this->seedInitialPotentials = (std::stoi(parameterValues.at(0)));
      // } else if(parameterName.find("seedInitialPrevSpike") != std::string::npos){
      //     userSeeds=true;
      //     this->seedInitialPreviousSpike = (std::stoi(parameterValues.at(0)));
    } else if (parameterName.find("noNeurons") != std::string::npos) {  //|| parameterName.find("nos") != std::string::npos
      this->noNeurons = std::stol(parameterValues.at(0));
    }
  }
  if (resetV >= thresholdV) {
    throw "There was an attempt to select a resetV>thresholdV. For safety reasons this is disabled.";
  }
  SetNeurons();
}

void NeuronPop::LoadPlasticityModel(const std::vector<FileEntry> &morphologyParameters) {
  if (!hasPlasticity) {
    for (auto &[parameterName, parameterValues] : morphologyParameters) {
      if (parameterName.find("pmodel_type") != std::string::npos) {
        if (parameterValues.at(0) == IDstringMonoDendriteSTDPTazerart) {
          morphologyType = IDstringMonoDendriteSTDPTazerart;
          for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
            (void)neuron;  // Does nothing, removes warning on unused vars
            this->morphology.push_back(std::make_unique<MonoDendriteSTDPTazerart>(
              this->infoGlobal,
              morphologyParameters));  // this conversion works but I do not remember why. Implicit downcasting through the move operation?
            this->hasPlasticity = true;
          }
        } else if (parameterValues.at(0) == IDstringMonoDendriteSTDPBiWindow) {
          morphologyType = IDstringMonoDendriteSTDPBiWindow;
          for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
            (void)neuron;  // Does nothing, removes warning on unused vars
            this->morphology.push_back(std::make_unique<MonoDendriteSTDPBiWindow>(this->infoGlobal, morphologyParameters));
            this->hasPlasticity = true;
            this->morphology.back()->LoadParameters(morphologyParameters);
          }
        } else if (parameterValues.at(0) == IDstringMonoDendriteSTDPTazerartRelative) {
          morphologyType = IDstringMonoDendriteSTDPTazerartRelative;
          for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
            (void)neuron;  // Does nothing, removes warning on unused vars
            this->morphology.push_back(std::make_unique<MonoDendriteSTDPTazerartRelative>(this->infoGlobal, morphologyParameters));
            this->hasPlasticity = true;
            this->morphology.back()->LoadParameters(morphologyParameters);
          }
        } else if (parameterValues.at(0) == IDstringTraceResourceHSTDP) {
          morphologyType = IDstringTraceResourceHSTDP;
          for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
            (void)neuron;  // Does nothing, removes warning on unused vars
            this->morphology.push_back(std::make_unique<AlphaResourceHSTDP>(this->infoGlobal, morphologyParameters));
            this->hasPlasticity = true;
            this->isBranched    = true;
            this->morphology.back()->LoadParameters(morphologyParameters);
          }

        } else if (parameterValues.at(0) == IDstringHeteroGraupnerBrunel) {
          morphologyType = IDstringHeteroGraupnerBrunel;
          for (NeuronInt neuron : std::ranges::views::iota(0, noNeurons)) {
            (void)neuron;  // Does nothing, removes warning on unused vars
            this->morphology.push_back(std::make_unique<HeteroGraupnerBrunel>(this->infoGlobal, morphologyParameters));
            this->hasPlasticity = true;
            this->isBranched    = true;
            this->morphology.back()->LoadParameters(morphologyParameters);
          }
        } else {
          throw "Unespecified morphology type";
        }
      }
    }
  } else {
    morphology.back()->CheckParameters(morphologyParameters);
  }
  if (!hasPlasticity) {  // This should never happen logically
    throw "Error, cannot use PModelSynapse without a plasticity model";
  }
}

void NeuronPop::SaveParameters(std::ofstream &wParameterStream) const {
  std::string neuronID = "neurons_" + std::to_string(GetId());

  wParameterStream << "#***********************************************\n";
  //*stream <<  neuronID + "_streamOutput                " << std::boolalpha<< streamingNOutputBool << std::noboolalpha << "\n";
  wParameterStream << neuronID + "_noNeurons\t\t\t" << std::to_string(noNeurons) << "\n";
  wParameterStream << neuronID + "_type\t\t\t\t" << GetType() << "\n";
  wParameterStream << neuronID + "_tauM\t\t\t\t" << std::to_string(this->membraneVTau) << " #secs\n";
  wParameterStream << neuronID + "_vReset\t\t\t" << std::to_string(this->resetV) << " #mV \n";
  wParameterStream << neuronID + "_vThresh\t\t\t" << std::to_string(this->thresholdV) << " #mV\n";
  wParameterStream << neuronID + "_refractoryTime\t\t" << std::to_string(this->refractorySteps * infoGlobal->dtTimestep) << " #secs\n";
  //*stream <<  neuronID + "_r_target                    " << std::to_string(this->targetRate)  << " Hz\n";
  if (userSeeds) {
    // wParameterStream <<  neuronID + "_seedInitialPotentials   " << this->seedInitialPotentials << "\n";
    // wParameterStream <<  neuronID + "_seedInitialPrevSpike    " << this->seedInitialPreviousSpike << "\n";
    wParameterStream << neuronID + "_seed\t\t\t\t" << this->seed << "\n";
  }
  if (definedPrevMeanFreq) {
    wParameterStream
      << neuronID + "_prevMinFrequency\t\t" << std::to_string(this->prevMinFrequency)
      << " #Hz #Average frequency at which the population was firing before the start of the simulation (relevant to STP classes).\n";
  }
  wParameterStream << "#\t\tNote: Resting potential is 0 by definition.\n";
}

void NeuronPop::SavePlasticityModel(std::ofstream &wParameterStream, std::string idString) const {
  if (this->hasPlasticity && (!this->savedModel)) {
    morphology.at(0)->SaveParameters(wParameterStream, idString);  //"neurons_" + std::to_string(GetId())
    return;
  } else if (savedModel) {
    return;
  } else {
    throw "This NeuronPop has no plasticity model";
  }
}

bool NeuronPop::HasSteadyState() const {
  if (morphology.empty()) {
    return false;
  } else {
    return morphology.at(0)->HasSteadyState();
  }
}

bool NeuronPop::ignoreJDParameters() const {
#ifndef NDEBUG
  if (morphology.empty()) {
    throw "Logical error jdis";
  } else {
    return morphology.at(0)->IgnoreJDParameters();
  }
#else
  return morphology.at(0)->IgnoreJDParameters();
#endif
}

// From here on

BaseSpinePtr NeuronPop::AllocateNewSynapse(NeuronInt neuronId, BranchTargeting &branchTarget) {
  std::lock_guard<std::mutex> _guardedMutexLock(_connectMutex);
  return morphology.at(neuronId)->AllocateNewSynapse(branchTarget);
}

std::string NeuronPop::GetIndividualSynapticProfileHeaderInfo() const {
  return morphology.at(0)->GetIndividualSynapticProfileHeaderInfo();
}

std::string NeuronPop::GetOverallSynapticProfileHeaderInfo() const {
  return morphology.at(0)->GetOverallSynapticProfileHeaderInfo();
}

std::vector<double> NeuronPop::GetIndividualSynapticProfile(NeuronInt neuronId, NeuronInt spineID) const {
  return this->morphology.at(neuronId)->GetIndividualSynapticProfile(spineID);
}

std::vector<double> NeuronPop::GetOverallSynapticProfile(NeuronInt neuronId) const {
  return this->morphology.at(neuronId)->GetOverallSynapticProfile();
}

std::vector<std::pair<std::string, double>> NeuronPop::GetSteadyStateData() const {
  return morphology.at(0)->GetSteadyStateData();
}

void NeuronPop::PostConnectSetUp() {
  std::lock_guard<std::mutex> _guardedMutexLock(_connectMutex);
  if (!morphology.empty()) {
    std::for_each(morphology.begin(), morphology.end(),
                  [](std::unique_ptr<Morphology> &singleMorphology) { singleMorphology->PostConnectSetUp(); });
  }
}