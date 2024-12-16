#include "NeuralNetwork.hpp"

NeuralNetwork::NeuralNetwork(std::string baseDirectory, std::vector<FileEntry> parameterEntries) {
  neurons  = nullptr;
  synapses = nullptr;
  recorder = nullptr;
  stimulus = nullptr;

  LoadParameters(baseDirectory, parameterEntries);

  SaveParameterOptions();
}

void NeuralNetwork::SaveParameters() const {
  // if this test file does not appear in the target directory: stop the
  // simulation and check the directoryPath.

  std::ofstream wParameterStream(recorder->GetParametersFilename());
  recorder->WriteHeader(wParameterStream);
  // stream <<  "#*****************************************************************\n";
  // stream <<  "#Date and Time:             " << std::ctime(&end_time);
  // stream <<  "#*****************************************************************\n";
  wParameterStream << "Title                       " << this->recorder->GetTitle() << "\n";
  wParameterStream << "#*****************************************************************\n";
  wParameterStream << "simulationTime              " << std::to_string(infoGlobal.simulationTime) << " \t\t#secs\n";
  wParameterStream << "dt_timestep                 " << std::to_string(infoGlobal.dtTimestep) << " \t\t#secs\n";
  wParameterStream << "globalSeed                  " << std::to_string(infoGlobal.globalSeed)
                   << " \t\t\t\t#overrides all other seeds if unequal -1\n";

  wParameterStream << "#*****************************************************************\n";
  wParameterStream << "#****************        Spatial parameters        ***************\n";
  wParameterStream << "#*****************************************************************\n";
  wParameterStream << "density                     " << std::to_string(infoGlobal.density)
                   << " \t\t\t\t#total number of neurons/mm^2 or /mm depending on Dimensions \n";
  wParameterStream << "Dimensions                  " << std::to_string(infoGlobal.dimensions) << " \t\t\t\t#for 1D set 1; for 2D set 2 \n";

  wParameterStream << "#*****************************************************************\n";
  wParameterStream << "#********** Scaling of synaptic and stimulus strengths ***********\n";
  wParameterStream << "#*****************************************************************\n";
  wParameterStream << "scalingSynapticStrength     " << std::to_string(infoGlobal.networkScaling_synStrength)
                   << "\t\t#Scaling exponent. Set = 0 if no scaling needed, otherwise typical exponent is -0.5.\n";
  wParameterStream << "scaling_C_N                 " << std::to_string(infoGlobal.networkScaling_mode)
                   << "\t\t\t\t# Set = 0 to scale with number of presynaptic neurons C. Set = 1 to scale with total number of neurons N. "
                      "(details below) \n";
  wParameterStream << "#\t\tscaling_C_N=0    scales internal synaptic strengths and UncorrelatedStimulus with C^s    and    "
                      "WhiteNoiseStimulus and SpatialGaussianStimulus with 1 \n";
  wParameterStream << "#\t\tscaling_C_N=1    scales internal synaptic strengths and UncorrelatedStimulus with N^s    and    "
                      "WhiteNoiseRescaled and SpatialGaussianStimulus with N^(-s) \n";
  wParameterStream
    << "#\t\tscalingSynapticStrength = s, N = number of neurons from all populations, C = average number of presynaptic neurons.\n";

  this->neurons->SaveParameters(wParameterStream);
  this->stimulus->SaveParameters(wParameterStream);
  this->recorder->SaveParameters(wParameterStream);
  this->synapses->SaveParameters(wParameterStream);

  wParameterStream << std::endl;
  wParameterStream.close();
}

void NeuralNetwork::SaveParameterOptions()
  const {  // This function should have stuff moved to the heap (and deleted at the end of the function)
  // return;//SaveParameteroptions is disabled
  std::cout << recorder->GetParameterOptionsFilename() << std::endl;
  std::ofstream streamPOptions(recorder->GetParameterOptionsFilename());

  GlobalSimInfo mockInfo      = infoGlobal;
  mockInfo.isMock             = true;
  GlobalSimInfo *mockInfo_ptr = &mockInfo;
  mockInfo_ptr->globalSeed =
    1;  // because many of the classes do not allow a negative seed, the ParameterOptions file is generated with an adapted globalInfo

  streamPOptions << "#************************************************************************************************\n";
  streamPOptions << "#*************  Stimulus options     ************************************************************\n";
  streamPOptions << "#************************************************************************************************\n";

  std::vector<FileEntry> stimulusString;
  std::string            placeholderString;

  placeholderString = "stimulus_meanCurrent ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "5.0 ";
  }
  stimulusString.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  placeholderString = "stimulus_sigmaCurrent ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "1.0 ";
  }
  stimulusString.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  std::unique_ptr<WhiteNoiseStimulus> whiteNoiseStimulus = std::make_unique<WhiteNoiseStimulus>(neurons, stimulusString, mockInfo_ptr);
  whiteNoiseStimulus->SaveParameters(streamPOptions);

  // std::unique_ptr<WhiteNoiseRescaled> whiteNoiseRescaledStimulus = std::make_unique<WhiteNoiseRescaled>(neurons,
  // stimulusString,mockInfo_ptr); whiteNoiseRescaledStimulus->SaveParameters(streamPOptions);

  placeholderString = "UncorrelatedStimulus_stimulus_J_X ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += std::to_string((static_cast<double>(neurons->GetTotalPopulations()) + 1) / 10) + " ";
  }
  stimulusString.push_back(SplitStringToEntry(placeholderString));
  stimulusString.push_back(
    SplitStringToEntry("UncorrelatedStimulus_stimulus_step " + std::to_string(infoGlobal.simulationTime) + " 5.00000"));
  std::unique_ptr<UncorrelatedPoissonLikeStimulus> uncorrelatedPoissonStimulus =
    std::make_unique<UncorrelatedPoissonLikeStimulus>(neurons, stimulusString, mockInfo_ptr);
  uncorrelatedPoissonStimulus->SaveParameters(streamPOptions);

  std::vector<FileEntry> whiteNoiseLinearStrings;
  placeholderString = "stimulus_meanCurrent ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "0.0 ";
  }
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "5.0 ";
  }
  whiteNoiseLinearStrings.push_back(SplitStringToEntry(placeholderString + "0.0 " + std::to_string(infoGlobal.simulationTime)));
  placeholderString = "stimulus_sigmaCurrent ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "0.0 ";
  }
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "1.0 ";
  }
  whiteNoiseLinearStrings.push_back(SplitStringToEntry(placeholderString + "0.0 " + std::to_string(infoGlobal.simulationTime)));
  std::unique_ptr<WhiteNoiseLinear> whiteNoiseLinearStimulus =
    std::make_unique<WhiteNoiseLinear>(neurons, whiteNoiseLinearStrings, mockInfo_ptr);
  whiteNoiseLinearStimulus->SaveParameters(streamPOptions);

  std::vector<FileEntry> spatialGaussianStimulusStrings;
  std::string            placeHolder;
  placeHolder = "stimulus_NumberOfGaussians           2";
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeHolder));
  placeHolder = "stimulus_X_position                  0.250000	 0.750000";
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeHolder));
  placeHolder = "stimulus_Y_position                  0.250000	 0.750000";
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeHolder));
  placeholderString = "stimulus_maxCurrent_1 ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "5.0 ";
  }
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  placeholderString = "stimulus_sigmaCurrent_t_1 ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "0.01 ";
  }
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  spatialGaussianStimulusStrings.push_back(
    SplitStringToEntry("stimulus_sigmaCurrent_x_1          1.0 " + std::to_string(infoGlobal.simulationTime)));
  placeholderString = "stimulus_maxCurrent_2 ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "5.0 ";
  }
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  placeholderString = "stimulus_sigmaCurrent_t_2 ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "0.01 ";
  }
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  spatialGaussianStimulusStrings.push_back(
    SplitStringToEntry("stimulus_sigmaCurrent_x_2          1.0 " + std::to_string(infoGlobal.simulationTime)));
  placeholderString = "stimulus_Background_Noise ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "1.00 ";
  }
  spatialGaussianStimulusStrings.push_back(SplitStringToEntry(placeholderString + std::to_string(infoGlobal.simulationTime)));
  std::unique_ptr<SpatialGaussianStimulus> SpatialGaussianStim =
    std::make_unique<SpatialGaussianStimulus>(neurons, spatialGaussianStimulusStrings, mockInfo_ptr);
  SpatialGaussianStim->SaveParameters(streamPOptions);

  stimulusString.push_back(SplitStringToEntry("stimulus_noExternalNeurons   1"));
  stimulusString.push_back(SplitStringToEntry("stimulus_PoissonTableEntries 10"));
  placeholderString = "ExternalConnectivity_lengthscale ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "0.25 ";
  }
  stimulusString.push_back(SplitStringToEntry(placeholderString));
  placeholderString = "ExternalConnectivity_PeakProbability ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    (void)neuronPop;
    placeholderString += "0.8 ";
  }
  stimulusString.push_back(SplitStringToEntry(placeholderString));
  placeholderString = "stimulus_J_X ";
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    placeholderString = placeholderString + std::to_string((static_cast<double>(neuronPop) + 1) / 10) + " ";
  }
  stimulusString.push_back(SplitStringToEntry(placeholderString));
  std::unique_ptr<SpatialPoissonStimulus> spatialPoissonStimulus =
    std::make_unique<SpatialPoissonStimulus>(neurons, stimulusString, mockInfo_ptr);
  spatialPoissonStimulus->SaveParameters(streamPOptions);

  streamPOptions << "#************************************************************************************************\n";
  streamPOptions << "#*************  Neuron options *****************************************************************\n";
  streamPOptions << "#************************************************************************************************\n";

  std::shared_ptr<LIFNeuronPop> neuronLIF = std::make_shared<LIFNeuronPop>(mockInfo_ptr, 0);
  neuronLIF->SaveParameters(streamPOptions);

  std::unique_ptr<EIFNeuronPop> neuronEIF = std::make_unique<EIFNeuronPop>(mockInfo_ptr, 0);
  neuronEIF->SaveParameters(streamPOptions);

  std::unique_ptr<QIFNeuronPop> neuronQIF = std::make_unique<QIFNeuronPop>(mockInfo_ptr, 0);
  neuronQIF->SaveParameters(streamPOptions);

  std::unique_ptr<PoissonNeuronPop> neuronPoisson = std::make_unique<PoissonNeuronPop>(mockInfo_ptr, 0);
  neuronPoisson->SaveParameters(streamPOptions);

  std::unique_ptr<DictatNeuronPop> neuronInput = std::make_unique<DictatNeuronPop>(mockInfo_ptr, 0);
  neuronInput->SaveParameters(streamPOptions);

  streamPOptions << "#************************************************************************************************\n";
  streamPOptions << "#*************  Connectivity options ************************************************************\n";
  streamPOptions << "#************************************************************************************************\n";

  std::string synapseString;
  synapseString                                   = "synapses_" + std::to_string(0) + "_" + std::to_string(0) + "_";
  std::shared_ptr<CurrentSynapse> synapse_current = std::make_shared<CurrentSynapse>(neuronLIF, neuronLIF, mockInfo_ptr);

  std::unique_ptr<RandomConnectivity> geometryRandomConnectivity =
    std::make_unique<RandomConnectivity>(static_cast<Synapse *>(synapse_current.get()), mockInfo_ptr);
  geometryRandomConnectivity->SaveParameters(streamPOptions, synapseString);

  streamPOptions << "#************************************************\n";
  std::unique_ptr<AdjacencyMatrixConnectivity> geometryAdjacencyMatrix =
    std::make_unique<AdjacencyMatrixConnectivity>(static_cast<Synapse *>(synapse_current.get()), mockInfo_ptr);
  geometryAdjacencyMatrix->SaveParameters(streamPOptions, synapseString);

  streamPOptions << "#************************************************\n";
  std::unique_ptr<PoissonConnectivity> geometryPoissonConnectivity =
    std::make_unique<PoissonConnectivity>(static_cast<Synapse *>(synapse_current.get()), mockInfo_ptr);
  geometryPoissonConnectivity->SaveParameters(streamPOptions, synapseString);
  streamPOptions << "#************************************************\n";
  std::unique_ptr<DistanceConnectivity> geometryDistanceConnectivity =
    std::make_unique<DistanceConnectivity>(static_cast<Synapse *>(synapse_current.get()), mockInfo_ptr);
  geometryDistanceConnectivity->SaveParameters(streamPOptions, synapseString);

  streamPOptions << "#************************************************************************************************\n";
  streamPOptions << "#*************  Synapse options *****************************************************************\n";
  streamPOptions << "#************************************************************************************************\n";

  synapse_current->SaveParameters(streamPOptions, synapseString);

  streamPOptions << "#************************************************\n";
  std::unique_ptr<MongilloSynapse> mongilloSynapse = std::make_unique<MongilloSynapse>(neuronLIF, neuronLIF, mockInfo_ptr);
  mongilloSynapse->SaveParameters(streamPOptions, synapseString);
  streamPOptions << "#************************************************\n";
  std::unique_ptr<MongilloSynapseContinuous> mongilloSynapseContinuous =
    std::make_unique<MongilloSynapseContinuous>(neuronLIF, neuronLIF, mockInfo_ptr);
  mongilloSynapseContinuous->SaveParameters(streamPOptions, synapseString);
  streamPOptions << "#************************************************\n";
  std::unique_ptr<PRGSynapseContinuous> PRGSContinuous = std::make_unique<PRGSynapseContinuous>(neuronLIF, neuronLIF, mockInfo_ptr);
  PRGSContinuous->SaveParameters(streamPOptions, synapseString);
  streamPOptions << "#************************************************\n";
  std::unique_ptr<PowerLawSynapse> powerLawSynapse = std::make_unique<PowerLawSynapse>(neuronLIF, neuronLIF, mockInfo_ptr);
  powerLawSynapse->SaveParameters(streamPOptions, synapseString);
  streamPOptions << "#************************************************\n";
  std::unique_ptr<ExponentialCurrentSynapse> exponentialCurrentSynapse =
    std::make_unique<ExponentialCurrentSynapse>(neuronLIF, neuronLIF, mockInfo_ptr);
  exponentialCurrentSynapse->SaveParameters(streamPOptions, synapseString);
  streamPOptions << "#************************************************\n";
  streamPOptions << "#************************************************************************************************\n";
  streamPOptions << "#*************  Heterosynaptic plasticity models ************************************************\n";
  streamPOptions << "#************************************************************************************************\n";
  streamPOptions << "#************************************************\n";
  std::vector<FileEntry>                    mockParams;
  std::unique_ptr<MonoDendriteSTDPTazerart> tazden = std::make_unique<MonoDendriteSTDPTazerart>(mockInfo_ptr, mockParams);
  tazden->SaveParameters(streamPOptions, "neurons_0");

  streamPOptions << "#*************  Branched dendrites **************************************************************\n";
  std::unique_ptr<AlphaResourceHSTDP> brhSTDP = std::make_unique<AlphaResourceHSTDP>(mockInfo_ptr, mockParams);
  brhSTDP->SaveParameters(streamPOptions, "neurons_0");
  // Instead of LIF, heteroLIF, put morphology options from Tazerart and HCS and HCP

  streamPOptions << std::endl;
  streamPOptions.close();
}

void NeuralNetwork::LoadParameters(std::string baseDirectory, std::vector<FileEntry> &parameterEntries) {
  // std::vector<std::string>        parameterEntry.parameterValues;
  // std::vector<std::string>        full_strs,neur_strs,syn_strs,stimulusString,rec_strs;
  std::string stimulusType;  //,neurons_type;
  std::string simulationTitle = "", nonIterateTitle;
  bool        timestepBool { false };

  for (FileEntry &parameterEntry : parameterEntries) {
    // parameterEntry.RemoveCommentsInValues();
    //  After this for loop, there is no need to check on # or other comments, they were forcefully removed in parameters
    if (parameterEntry.parameterName.compare(TitleIDString) == 0) {
      if (parameterEntry.parameterValues.empty()) {  // Should never happen
        simulationTitle = "default";
        // nonIterateTitle = "default";
      } else {
        simulationTitle = parameterEntry.parameterValues.at(0);
        // nonIterateTitle = parameterEntry.parameterValues.at(1);
      }
    } else if ((parameterEntry.parameterName.find("dt_timestep") != std::string::npos ||
                parameterEntry.parameterName.find("dt") != std::string::npos) &&
               !timestepBool) {
      infoGlobal.dtTimestep = std::stod(parameterEntry.parameterValues.at(0));
      infoGlobal.dtSqrt     = sqrt(infoGlobal.dtTimestep);
      timestepBool          = true;
    } else if (parameterEntry.parameterName.find("simulationTime") != std::string::npos) {
      infoGlobal.simulationTime = std::stod(parameterEntry.parameterValues.at(0));
    } else if (parameterEntry.parameterName.find("globalSeed") != std::string::npos) {
      infoGlobal.globalSeed      = std::stoi(parameterEntry.parameterValues.at(0));
      infoGlobal.globalGenerator = std::mt19937(infoGlobal.globalSeed);
    } else if (parameterEntry.parameterName.find("density") != std::string::npos) {
      infoGlobal.density = std::stoi(parameterEntry.parameterValues.at(0));
    } else if (parameterEntry.parameterName.find("Dimensions") != std::string::npos) {
      infoGlobal.dimensions = std::stoi(parameterEntry.parameterValues.at(0));
    } else if ((parameterEntry.parameterName.find("synapticScaling") != std::string::npos) ||
               (parameterEntry.parameterName.find("scalingSynapticStrength") != std::string::npos)) {
      infoGlobal.networkScaling_synStrength = std::stod(parameterEntry.parameterValues.at(0));
    } else if (parameterEntry.parameterName.find("scaling_C_N") != std::string::npos) {
      infoGlobal.networkScaling_mode = std::stoi(parameterEntry.parameterValues.at(0));
      // } else if(parameterEntry.parameterName.find("neurons_type") != std::string::npos){
      //     neurons_type = parameterEntry.parameterValues.at(0);//Is this still relevant ?
      // } else if(parameterEntry.parameterName.find("recorder_type") != std::string::npos){
      //     recorderType = parameterEntry.parameterValues.at(0);
    } else if (parameterEntry.parameterName.find("stimulus_type") != std::string::npos) {
      stimulusType = parameterEntry.parameterValues.at(0);
    } else if (parameterEntry.parameterName.find("pathToInputFile") != std::string::npos) {
      if (parameterEntry.parameterValues.size() == 0) {
        infoGlobal.pathToInputFile = "";
      } else {
        infoGlobal.pathToInputFile = parameterEntry.parameterValues.at(0);
      }
    } else if (parameterEntry.parameterName.find(NonIterateTitleIDString) != std::string::npos) {
      nonIterateTitle = parameterEntry.parameterValues.at(0);
    }
  }
  infoGlobal.pathToInputFile = infoGlobal.pathToInputFile + nonIterateTitle + "_";

  std::vector<FileEntry> neuronParameters { FilterStringEntries(parameterEntries, "neurons_") };
  std::vector<FileEntry> synapseParameters { FilterStringEntries(parameterEntries, "synapses_") };
  std::vector<FileEntry> stimulusParameters { FilterStringEntries(parameterEntries, "stimulus_") };
  std::vector<FileEntry> recorderParameters { FilterStringEntries(parameterEntries, "recorder_") };

  // assemble vector of lines for classes that are not yet adapted to ParameterEntry structs
  //  for (auto & parEntry : *parEntries) {
  //      line.clear();
  //      line = parEntry.parameterEntry.parameterName + " ";
  //      for(auto & value : parEntry.parameterEntry.parameterValues)
  //          line.append(value+" ");
  //      full_strs.push_back(line);
  //  }

  // FilterStringVector(&full_strs,"neurons_",  &neur_strs);
  // FilterStringVector(&full_strs,"synapses_", &syn_strs);
  // FilterStringVector(&full_strs,"stimulus_", &stimulusString);
  // FilterStringVector(&full_strs,"recorder_", &rec_strs);

  /*if(neurons_type == stringLIFNeuron)
      this->neurons   = new LIFNeuronPop(&infoGlobal);
  else if(neurons_type == stringEIFNeuron)
      this->neurons   = new EIFNeuronPop(&infoGlobal);
  else if(neurons_type == stringQIFNeuron)
      this->neurons   = new QIFNeuronPop(&infoGlobal);
  else{
      std::cout << "Neuron not defined";
      return 0;
   }*/

  neurons  = std::make_shared<NeuronPopSample>(neuronParameters, &infoGlobal);
  synapses = std::make_shared<SynapseSample>(neurons, synapseParameters, &infoGlobal);

  if (stimulusType == IDstringUncorrelatedStimulus) {
    this->stimulus = std::make_shared<UncorrelatedPoissonLikeStimulus>(neurons, stimulusParameters, &infoGlobal);  //
  } else if (stimulusType == IDstringWhitenoiseStimulus) {
    this->stimulus = std::make_shared<WhiteNoiseStimulus>(neurons, stimulusParameters, &infoGlobal);
  } else if (stimulusType == IDstringWhitenoiseRescaled) {
    this->stimulus = std::make_shared<WhiteNoiseStimulus>(neurons, stimulusParameters, &infoGlobal);  //
    if (this->infoGlobal.networkScaling_mode == 0) {
      throw "WhiteNoiseRescaled (legacy) cannot have scaling_C_N = 0";
    }
  } else if (stimulusType == IDstringWhiteNoiseLinear) {
    this->stimulus = std::make_shared<WhiteNoiseLinear>(neurons, stimulusParameters, &infoGlobal);  //
  } else if (stimulusType == IDstringSpatialGaussianStimulus) {
    this->stimulus = std::make_shared<SpatialGaussianStimulus>(neurons, stimulusParameters, &infoGlobal);  //
  } else if (stimulusType == IDstringSpatialPoissonStimulus) {
    this->stimulus = std::make_shared<SpatialPoissonStimulus>(neurons, stimulusParameters, &infoGlobal);  //
  } else if (stimulusType == "" || stimulusType == IDstringNoStimulus) {
    this->stimulus = std::make_shared<NoStimulus>(neurons, &infoGlobal);
  } else {
    throw "Stimulus type was not properly defined.";
  }

  this->recorder = std::make_unique<Recorder>(neurons, synapses, stimulus, baseDirectory, recorderParameters, simulationTitle,
                                              nonIterateTitle, &infoGlobal);
}

bool NeuralNetwork::WellDefined() const {
  return !((neurons == nullptr) || (synapses == nullptr) || (recorder == nullptr) || (stimulus == nullptr));
}

void NeuralNetwork::Simulate() {
  if (!WellDefined()) {
    throw "NeuralNetwork was not properly defined";
  }

  //******************************
  // Declarations & Initialization
  //******************************
  double                                 timestepProgressRatio { 0.0 };
  std::chrono::seconds                   computationTime;
  PopInt                                 totalNeuronPops { neurons->GetTotalPopulations() };
  std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
  TStepInt simSteps { static_cast<TStepInt>(infoGlobal.simulationTime / infoGlobal.dtTimestep) };  // number of simulation time steps
  // int      global_D_max = this->synapses->GetMaxD();          // get maximum delay across all synapses: size of waiting matrix DEPRECATED

  std::vector<std::vector<double>> synaptic_dV(totalNeuronPops);
  for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
    synaptic_dV.at(neuronPop).resize(neurons->GetNeuronsPop(neuronPop));
  }

  this->recorder->SetFilenameDate();
  SaveParameters();
  //*****************************************************
  // --------------- SIMULATION SETUP ------------
  //*****************************************************
  infoGlobal.timeStep = 0;
  // infoGlobal.waitingIndex = 0;
  auto setupStart { std::chrono::high_resolution_clock::now() };
  auto intermediateTime { setupStart };
  // Setup
  this->synapses->ConnectNeurons();
  this->recorder->WriteConnectivity();
  this->recorder->WriteDistributionD();
  this->recorder->WriteDistributionJ();
  this->recorder->WriteDataHeader();
  //*****************************************************
  // --------------- START OF THE SIMULATION ------------
  //*****************************************************
  std::cout << "\n (Non)Pandas start simulation : " << this->recorder->GetTitle() << std::endl;
  auto simulateStart = std::chrono::high_resolution_clock::now();
  // std::chrono::microseconds duration;
  // Remember remember
  // The 5th of November
  while (infoGlobal.timeStep <= simSteps) {
    // auto t1 = std::chrono::high_resolution_clock::now();
    for (std::vector<double> &popSynaptic_dV : synaptic_dV) {
      std::fill(popSynaptic_dV.begin(), popSynaptic_dV.end(), 0.0);
    }
    // auto t2 = std::chrono::high_resolution_clock::now();
    this->synapses->Advect(synaptic_dV);
    // auto t3 = std::chrono::high_resolution_clock::now();
    this->stimulus->Update(synaptic_dV);
    // auto t4 = std::chrono::high_resolution_clock::now();
    this->neurons->Advect(synaptic_dV);  // spikers renewed
    // auto t5 = std::chrono::high_resolution_clock::now();
    this->recorder->Record(synaptic_dV);
    // auto t6 = std::chrono::high_resolution_clock::now();
    this->synapses->Reset();

    // duration = duration_cast<std::chrono::microseconds>(t2 - t1);
    // std::cout<<"synaptic_dV clearing:" << duration.count() << std::endl;
    // duration = duration_cast<std::chrono::microseconds>(t3 - t2);
    // std::cout<< "synapse advect:" << duration.count()<<std::endl;
    // duration = duration_cast<std::chrono::microseconds>(t4 - t3);
    // std::cout<< "Stimulus update:" << duration.count()<<std::endl;
    // duration = duration_cast<std::chrono::microseconds>(t5 - t4);
    // std::cout<< "Neuron advect:" << duration.count()<<std::endl;
    // duration = duration_cast<std::chrono::microseconds>(t6 - t5);
    // std::cout<< "Recorder Record:" << duration.count()<<std::endl;
    /*        double accumulator{};
    for (std::vector<double>& popSynaptic_dV : synaptic_dV) {
        accumulator = std::accumulate(popSynaptic_dV.begin(), popSynaptic_dV.end(), accumulator, [](double accumulator, double value)
    {return accumulator + value; });
    }
    std::cout << "accumulated syndv: "<<accumulator <<"\n";*/

    if (infoGlobal.timeStep % (static_cast<long>(simSteps * 0.001)) == 1) {
      intermediateTime      = std::chrono::high_resolution_clock::now();
      timestepProgressRatio = (static_cast<double>(infoGlobal.timeStep) / static_cast<double>(simSteps));
      computationTime       = duration_cast<std::chrono::seconds>(intermediateTime - setupStart);
      std::cout << std::fixed << std::setprecision(1) << timestepProgressRatio * 100 << "%  -- Comp. time: " << computationTime.count()
                << "/" << static_cast<int>(computationTime.count() / timestepProgressRatio) << " sec. -- " << std::endl;
    }

    infoGlobal.timeStep++;
  }

  this->recorder->Record(synaptic_dV);
  auto simulateEnd    = std::chrono::high_resolution_clock::now();
  auto simulationTime = std::chrono::duration_cast<std::chrono::seconds>(simulateEnd - simulateStart);
  auto setupTime      = std::chrono::duration_cast<std::chrono::seconds>(simulateStart - setupStart);
  this->recorder->WriteFinalDataFile(setupTime, simulationTime);
  std::cout << "\n(Not)Pandas end simulation : " << this->recorder->GetTitle() << "\n";
  std::cout << "\nTotal simulation time(s): " << simulationTime.count() << "\n";

  // We should also make sure this content is stored in an output file .txt. This function is a console output for debugging.

  //*****************************************************
  // --------------- END OF THE SIMULATION ------------
  //*****************************************************
}

void NeuralNetwork::MakeInputCopies(const std::string &inputFileAddress) const {
  this->recorder->MakeInputCopies(inputFileAddress);
}