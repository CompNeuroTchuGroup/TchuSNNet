
#include "WhiteNoiseLinear.hpp"

WhiteNoiseLinear::WhiteNoiseLinear(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &stimulusParameters, GlobalSimInfo *infoGlobal)
    : Stimulus(neurons, infoGlobal) {

  //******************************
  //***** Default parameterValues *********
  //******************************
  //******************************
  cachedScalingConstants.resize(neurons->GetTotalPopulations());
  meanScalingCache.resize(neurons->GetTotalPopulations());
  sigmaScalingCache.resize(neurons->GetTotalPopulations());
  LoadParameters(stimulusParameters);
  // generator = std::mt19937(seed);
}

void WhiteNoiseLinear::LoadParameters(const std::vector<FileEntry> &stimulusParameters) {

  Stimulus::LoadParameters(stimulusParameters);
  PopInt totalNeuronPops = neurons->GetTotalPopulations();

  for (auto &[parameterName, parameterValues] : stimulusParameters) {
    // if ((parameterName.find("seed") != std::string::npos)) {
    // 	seed = static_cast<signed int>(std::stod(parameterValues.at(0)));
    // }
    if ((parameterName.find("meanCurrent") != std::string::npos)) {
      StepStruct step;
      for (PopInt readingIndex : std::ranges::views::iota(0, totalNeuronPops)) {
        // Here we are only reading the input parameters, which they should be 2*(totalNeuronPops+1). Their order of usage will differ of order of
        // reading
        step.parameterValues.push_back(std::stod(parameterValues.at(2 * static_cast<size_t>(readingIndex))));
        step.parameterValues.push_back(std::stod(parameterValues.at(2 * static_cast<size_t>(readingIndex) + 1)));
      }
      step.startTimeStep = static_cast<int>(std::round(std::stod(parameterValues.at(2 * totalNeuronPops)) / infoGlobal->dtTimestep));
      step.endTimeStep   = static_cast<int>(std::round(std::stod(parameterValues.at(2 * totalNeuronPops + 1)) / infoGlobal->dtTimestep));

      meanCurrent.push_back(step);
    } else if ((parameterName.find("sigmaCurrent") != std::string::npos)) {
      StepStruct step;
      for (PopInt readingIndex : std::ranges::views::iota(0, totalNeuronPops)) {
        step.parameterValues.push_back(std::stod(parameterValues.at(2 * readingIndex)));
        step.parameterValues.push_back(std::stod(parameterValues.at(2 * readingIndex + 1)));
      }
      step.startTimeStep = static_cast<int>(std::round(std::stod(parameterValues.at(2 * totalNeuronPops)) / infoGlobal->dtTimestep));
      step.endTimeStep   = static_cast<int>(std::round(std::stod(parameterValues.at(2 * totalNeuronPops + 1)) / infoGlobal->dtTimestep));
      sigmaCurrent.push_back(step);
    }
  }
  if (infoGlobal->isMock) {
    return;
  }
  PostLoadParameters();
}

void WhiteNoiseLinear::SaveParameters(std::ofstream &wParameterStream) const {

  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  Stimulus::SaveParameters(wParameterStream);

  // if (infoGlobal->globalSeed == -1) {
  // 	wParameterStream << "stimulus_seed                        " << std::to_string(seed) << "\n";
  // }

  for (const StepStruct &step : meanCurrent) {
    wParameterStream << "stimulus_meanCurrent                 ";
    for (PopInt writingIndex : std::ranges::views::iota(0, totalNeuronPops)) {
      wParameterStream << std::to_string(step.parameterValues.at(2 * writingIndex)) << "\t ";
      wParameterStream << std::to_string(step.parameterValues.at(2 * writingIndex + 1)) << "\t ";
    }
    wParameterStream << std::to_string(static_cast<double>(step.startTimeStep) * infoGlobal->dtTimestep) << " \t";
    wParameterStream << std::to_string(static_cast<double>(step.endTimeStep) * infoGlobal->dtTimestep) << " \t";
    wParameterStream
        << " #[col 1: input to pop0 at t_0, col 2: pop1 at t_0, ... colP+1: pop1 t_f, ... col2P: popN t_f, t0, tf. Dimensions: [mV/sec , secs.]\n";
  }

  for (const StepStruct &step : sigmaCurrent) {
    wParameterStream << "stimulus_sigmaCurrent                ";
    for (PopInt writingIndex : std::ranges::views::iota(0, totalNeuronPops)) {
      wParameterStream << std::to_string(step.parameterValues.at(2 * writingIndex)) << "\t ";
      wParameterStream << std::to_string(step.parameterValues.at(2 * writingIndex + 1)) << "\t ";
    }
    wParameterStream << std::to_string(static_cast<double>(step.startTimeStep) * infoGlobal->dtTimestep) << " \t";
    wParameterStream << std::to_string(static_cast<double>(step.endTimeStep) * infoGlobal->dtTimestep) << " \t";
    wParameterStream
        << " #[col 1: input to pop0 at t_0, col 2: pop1 at t_0, ... colP+1: pop1 t_f, ... col2P: popN t_f, t0, tf. Dimensions: [mV/sqrt(sec) , secs.]\n";
  }

  wParameterStream << "#\t\tRI_{i,ext}/tauM*dt = meanCurrent_i*dt + sqrt(dt)*sigmaCurrent_i*NormalDistribution(0,1)\n";
}

void WhiteNoiseLinear::PostLoadParameters() {
  std::cout << "\nSetting up the stimulus class...";
  // if(meanCurrent.empty()){ //Nothing should happen if there is no structs in the vector
  //     return;
  // }
  if (meanCurrent.empty() || sigmaCurrent.empty()) {
    wellDefined = false;
    std::cout << "Stimulus class was ill-defined\n" << std::endl;
    return;
  }
  PopInt totalPopulations{neurons->GetTotalPopulations()};
  // Should we sort the vectors here? This assumes they are in order
  meanCurrentPtr  = meanCurrent.data();
  sigmaCurrentPtr = sigmaCurrent.data();
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    cachedScalingConstants.at(neuronPop) = infoGlobal->dtTimestep * pow(GetScaling(neuronPop), -(infoGlobal->networkScaling_synStrength));
  }
  for (PopInt neuronPop : std::ranges::views::iota(0, totalPopulations)) {
    meanScalingCache.at(neuronPop) =
        (meanCurrentPtr->parameterValues.at(neuronPop + totalPopulations) - meanCurrentPtr->parameterValues.at(neuronPop)) /
        (meanCurrentPtr->endTimeStep - meanCurrentPtr->startTimeStep);
  }
  for (PopInt neuronPop : std::ranges::views::iota(0, totalPopulations)) {
    sigmaScalingCache.at(neuronPop) =
        (sigmaCurrentPtr->parameterValues.at(neuronPop + totalPopulations) - sigmaCurrentPtr->parameterValues.at(neuronPop)) /
        (sigmaCurrentPtr->endTimeStep - sigmaCurrentPtr->startTimeStep);
  }
  std::cout << "Done.\n" << std::endl;
}

void WhiteNoiseLinear::SetSignalMatrix() {

  // Refactor whole function

  double mean, sigma;
  double initialTime;
  double initialMean, initialSigma;
  // double dtTimestep = infoGlobal->dtTimestep;
  TStepInt timeStep = infoGlobal->timeStep;

  double dtSqrtd{infoGlobal->dtSqrt};
  // Pointer arithmetic is unsafe, find an alternative
  //  StepStruct   *mean_step_current = &meanCurrent.at(0);
  //  StepStruct   *stepCurrentSigma = &sigmaCurrent.at(0);

  // Rewrite these two while loops with std::find_if, making sure the last iterator given derefernces to the last instruction.
  if ((meanCurrentPtr != &meanCurrent.back()) && (meanCurrentPtr->endTimeStep <= timeStep)) {
    meanCurrentPtr++;
    RecalculateMeanScalingCache(meanCurrentPtr);
  }
  // const StepStruct& meanStepCurrent = *std::find_if(meanCurrent.begin(), meanCurrent.end(), [timeStep](StepStruct step){
  // //Same as the while statement
  // 		return (timeStep<=step.endTimeStep) || (step.lastStep);
  // });
  if ((sigmaCurrentPtr != &sigmaCurrent.back()) && (sigmaCurrentPtr->endTimeStep <= timeStep)) {
    sigmaCurrentPtr++;
    RecalculateSigmaScalingCache(sigmaCurrentPtr);
  }
  // const StepStruct& stepCurrentSigma = *std::find_if(sigmaCurrent.begin(), sigmaCurrent.end(), [timeStep](StepStruct step){
  // //Same as the while statement
  // 		return (timeStep<=step.endTimeStep) || (step.lastStep);
  // });
  // while ((timeStep > stepCurrentSigma->endTimeStep) && (stepCurrentSigma != &sigmaCurrent.back()))
  // 	stepCurrentSigma++;

  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    // Order of indexing is according to what is stated in the comment of the parameter
    initialMean = meanCurrentPtr->parameterValues.at(neuronPop);
    initialTime = meanCurrentPtr->startTimeStep;
    mean        = initialMean + meanScalingCache.at(neuronPop) * (timeStep - initialTime);
    // All of these calculations could have been made before the Simulate loop and stored inside the struct
    initialSigma = sigmaCurrentPtr->parameterValues.at(neuronPop);
    initialTime  = sigmaCurrentPtr->startTimeStep;
    sigma        = initialSigma + sigmaScalingCache.at(neuronPop) * (timeStep - initialTime);
    for (NeuronInt neuron : std::ranges::views::iota(0, neurons->GetNeuronsPop(neuronPop))) {
      signalMatrix.at(neuronPop).at(neuron) = mean * cachedScalingConstants.at(neuronPop) + dtSqrtd * sigma * standardDistribution(generator);
      // this is constant and independent of the rest of the simulation. This could be precalculated but it would require a lot of memory (one matrix
      // per timestep)
    }
  }
}

void WhiteNoiseLinear::RecalculateMeanScalingCache(StepStruct *meanStep) {
  PopInt totalPopulations{neurons->GetTotalPopulations()};
  for (PopInt neuronPop : std::ranges::views::iota(0, totalPopulations)) {
    meanScalingCache.at(neuronPop) = (meanStep->parameterValues.at(neuronPop + totalPopulations) - meanStep->parameterValues.at(neuronPop)) /
                                     (meanStep->endTimeStep - meanStep->startTimeStep);
  }
}

void WhiteNoiseLinear::RecalculateSigmaScalingCache(StepStruct *sigmaStep) {
  PopInt totalPopulations{neurons->GetTotalPopulations()};
  for (PopInt neuronPop : std::ranges::views::iota(0, totalPopulations)) {
    sigmaScalingCache.at(neuronPop) = (sigmaStep->parameterValues.at(neuronPop + totalPopulations) - sigmaStep->parameterValues.at(neuronPop)) /
                                      (sigmaStep->endTimeStep - sigmaStep->startTimeStep);
  }
}

void WhiteNoiseLinear::Update(std::vector<std::vector<double>> &synaptic_dV) {
  SetSignalMatrix();
  Stimulus::Update(synaptic_dV);
}

double WhiteNoiseLinear::GetScaling(PopInt neuronPop) const {
  if (infoGlobal->networkScaling_mode == 0)
    return 1.0;
  else if (infoGlobal->networkScaling_mode == 1)
    return (neurons->GetTotalNeurons());
  else {
    throw "WhiteNoiseLinear (GetScaling): External Scaling not well defined! (0 or 1)";
  }
};
