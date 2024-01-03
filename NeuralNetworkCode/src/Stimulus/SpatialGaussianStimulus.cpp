
//  Created by Pierre Ekelmans on 20/03/2018.
//  Copyright © 2018 Pierre Ekelmans. All rights reserved.
//

#include "SpatialGaussianStimulus.hpp"

SpatialGaussianStimulus::SpatialGaussianStimulus(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &stimulusParameters,
                                                 GlobalSimInfo *infoGlobal)
    : Stimulus(neurons, infoGlobal) {

  //******************************
  //***** Default parameterValues *********
  //******************************

  LoadParameters(stimulusParameters);
}

void SpatialGaussianStimulus::LoadParameters(const std::vector<FileEntry> &stimulusParameters) {

  Stimulus::LoadParameters(stimulusParameters);
  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  int    gaussianIndex;

  maxCurrent.resize(noGaussians);
  sigmaCurrentT.resize(noGaussians);
  sigmaCurrentX.resize(noGaussians);
  gaussianPositionX.resize(noGaussians, 0.5);
  if (infoGlobal->dimensions == 2) {
    gaussianPositionY.resize(noGaussians, 0.5);
  } else {
    gaussianPositionY.resize(noGaussians);
  }

  for (auto &[parameterName, parameterValues] : stimulusParameters) {
    // By default, unless n is defined (see following input load), there is only one Gaussian, centered in 0.5 0.5
    if ((parameterName.find("NumberOfGaussians") != std::string::npos)) {
      noGaussians = std::stoi(parameterValues.at(0));
      /*			if (Ngauss > 9) {
                      std::cout << "Only up to 9 Gaussians are supported at this point";
                      Ngauss = 9;
                  }*/
      maxCurrent.resize(noGaussians);
      sigmaCurrentT.resize(noGaussians);
      sigmaCurrentX.resize(noGaussians);
      maxCurrentPtrs.resize(noGaussians);
      sigmaCurrentTPtrs.resize(noGaussians);
      sigmaCurrentXPtrs.resize(noGaussians);
      gaussianPositionX.resize(noGaussians);
      gaussianPositionY.resize(noGaussians);
      cachedMuNeuronFactors.resize(noGaussians);
      for (std::vector<std::vector<double>> &gaussianCache : cachedMuNeuronFactors) {
        gaussianCache.resize(neurons->GetTotalPopulations());
        for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
          gaussianCache.at(neuronPop).resize(neurons->GetNeuronsPop(neuronPop));
        }
      }

    } else if ((parameterName.find("X_position") != std::string::npos)) {
      for (int gaussian : std::ranges::views::iota(0, noGaussians)) {
        gaussianPositionX.at(gaussian) = std::stod(parameterValues.at(gaussian));
      }
    } else if ((parameterName.find("Y_position") != std::string::npos) && (infoGlobal->dimensions == 2)) {
      for (int gaussian : std::ranges::views::iota(0, noGaussians)) {
        gaussianPositionY.at(gaussian) = std::stod(parameterValues.at(gaussian));
      }
    } else if ((parameterName.find("maxCurrent") != std::string::npos)) {
      StepStruct step;
      gaussianIndex = 1;
      if (parameterName.find("t_") != std::string::npos) {
        gaussianIndex = std::stoi(parameterName.substr(parameterName.find("t_") + 2, std::string::npos));
      }
      step.parameterValues.resize(totalNeuronPops);
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        step.parameterValues.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
      }
      if (isDouble(parameterValues.at(totalNeuronPops))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
      } else {
        step.endTimeStep = LONG_MAX;
      }
      if ((step.endTimeStep < 0) || (step.endTimeStep > static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep)))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep));
      }
      maxCurrent.at(gaussianIndex - 1).push_back(step);
    } else if ((parameterName.find("sigmaCurrent_x") != std::string::npos)) {
      StepStruct step;
      gaussianIndex = 1;
      if (parameterName.find("x_") != std::string::npos) {
        gaussianIndex = std::stoi(parameterName.substr(parameterName.find("x_") + 2, std::string::npos));
      }
      step.parameterValues.resize(1);
      step.parameterValues.at(0) = std::stod(parameterValues.at(0));
      if (isDouble(parameterValues.at(1))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(std::stod(parameterValues.at(1)) / infoGlobal->dtTimestep));
      } else {
        step.endTimeStep = LONG_MAX;
      }
      if ((step.endTimeStep < 0) || (step.endTimeStep > static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep)))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep));
      }
      sigmaCurrentX.at(gaussianIndex - 1).push_back(step);
    } else if ((parameterName.find("sigmaCurrent_t") != std::string::npos)) {
      StepStruct step;
      gaussianIndex = 1;
      if (parameterName.find("_t_") != std::string::npos) {
        gaussianIndex = std::stoi(parameterName.substr(parameterName.find("_t_") + 3, std::string::npos));
      }
      step.parameterValues.resize(totalNeuronPops);
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        step.parameterValues.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
      }
      if (isDouble(parameterValues.at(totalNeuronPops))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
      } else {
        step.endTimeStep = LONG_MAX;
      }
      if ((step.endTimeStep < 0) || (step.endTimeStep > static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep))))
        step.endTimeStep = static_cast<TStepInt>(std::round(infoGlobal->simulationTime / infoGlobal->dtTimestep));
      sigmaCurrentT.at(gaussianIndex - 1).push_back(step);
    } else if ((parameterName.find("Background_Noise") != std::string::npos)) {
      StepStruct step;
      step.parameterValues.resize(totalNeuronPops);
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        step.parameterValues.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
      }
      if (isDouble(parameterValues.at(totalNeuronPops))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
      } else {
        step.endTimeStep = LONG_MAX;
      }
      if ((step.endTimeStep < 0) || (step.endTimeStep > static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep)))) {
        step.endTimeStep = static_cast<TStepInt>(std::lround(infoGlobal->simulationTime / infoGlobal->dtTimestep));
      }
      backgroundNoise.push_back(step);
    }
  }
  if (infoGlobal->isMock) {
    return;
  }
  PostLoadParameters();
}

void SpatialGaussianStimulus::SaveParameters(std::ofstream &wParameterStream) const {

  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  Stimulus::SaveParameters(wParameterStream);

  wParameterStream << "stimulus_NumberOfGaussians           ";
  wParameterStream << std::to_string(noGaussians) << "\n";

  wParameterStream << "stimulus_X_position                  ";
  for (int gaussianIndex : std::ranges::views::iota(0, noGaussians)) {
    wParameterStream << std::to_string(gaussianPositionX.at(gaussianIndex)) << "\t ";
  }
  wParameterStream << "#Position of each Gaussian on the X axis (between 0 and 1)\n";

  if (infoGlobal->dimensions == 2) {
    wParameterStream << "stimulus_Y_position                  ";
    for (int gaussianIndex : std::ranges::views::iota(0, noGaussians)) {
      wParameterStream << std::to_string(gaussianPositionY.at(gaussianIndex)) << "\t ";
    }
    wParameterStream << "#Position of each Gaussian on the Y axis (between 0 and 1)\n";
  }

  for (int gaussianIndex : std::ranges::views::iota(0, noGaussians)) {
    for (const StepStruct &step : maxCurrent.at(gaussianIndex)) {
      wParameterStream << "stimulus_maxCurrent_" << std::to_string(gaussianIndex + 1) << "                ";
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        wParameterStream << std::to_string(step.parameterValues.at(neuronPop)) << "\t ";
      }
      wParameterStream << std::to_string(static_cast<double>(step.endTimeStep) * infoGlobal->dtTimestep) << " \t";
      wParameterStream
          << " #[column i: input to neurons of population i, at the center of the gaussian, last column: time until which input is set. Dimensions: [mV/sec , secs.]\n";
    }
    for (const StepStruct &step : sigmaCurrentT.at(gaussianIndex)) {
      wParameterStream << "stimulus_sigmaCurrent_t_" << std::to_string(gaussianIndex + 1) << "            ";
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        wParameterStream << std::to_string(step.parameterValues.at(neuronPop)) << "\t ";
      }
      wParameterStream << std::to_string(static_cast<double>(step.endTimeStep) * infoGlobal->dtTimestep) << " \t";
      wParameterStream
          << " #[column i: relative input noise to population i (relative to the mean current), last column: time until which input is set. Dimensions: [ -  , secs.]\n";
    }

    for (const StepStruct &step : sigmaCurrentX.at(gaussianIndex)) {
      wParameterStream << "stimulus_sigmaCurrent_x_" << std::to_string(gaussianIndex + 1) << "            ";
      wParameterStream << std::to_string(step.parameterValues.at(0)) << "\t "; // the width of the input is the same for all popluations
      wParameterStream << std::to_string(static_cast<double>(step.endTimeStep) * infoGlobal->dtTimestep) << " \t";
      wParameterStream
          << " #[column 1: spatial spread (std of the Gaussian) of the input to all populations, last column: time until which input is set. Dimensions: [mm , secs.]\n";
    }
  }

  for (const StepStruct &step : backgroundNoise) {
    wParameterStream << "stimulus_Background_Noise            ";
    for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
      wParameterStream << std::to_string(step.parameterValues.at(neuronPop)) << "\t ";
    }
    wParameterStream << std::to_string(static_cast<double>(step.endTimeStep) * infoGlobal->dtTimestep) << " \t";
    wParameterStream << "#Noise applied in the whole domain [mV/sqrt(sec) , secs.]\n";
  }
  wParameterStream
      << "#\t\tRI_{i,ext}/tauM*dt = meanCurrent_i*dt*exp(-d{i}^2/(2sigmaCurrent_x)) + sqrt(dt)*sigmaCurrent_t_i*NormalDistribution(0,1)\t Where d{i} is the distance of neuron i to the center of the domain\n";
}
void SpatialGaussianStimulus::PostLoadParameters() {
  std::cout << "\nSetting up the stimulus class...";
  // Well defined checks
  auto stepMatrixEmptyBool = [](std::vector<StepStruct> vector) { return vector.empty(); };
  if (backgroundNoise.empty() || std::any_of(sigmaCurrentX.begin(), sigmaCurrentX.end(), stepMatrixEmptyBool) ||
      std::any_of(sigmaCurrentT.begin(), sigmaCurrentT.end(), stepMatrixEmptyBool) ||
      std::any_of(maxCurrent.begin(), maxCurrent.end(), stepMatrixEmptyBool)) {
    wellDefined = false;
    std::cout << "Stimulus class was ill-defined\n" << std::endl;
    return;
  }
  // Set up pointer structures
  backgroundNoisePtr = backgroundNoise.data();
  for (int gaussianIndex : std::ranges::views::iota(0, noGaussians)) {
    maxCurrentPtrs.at(gaussianIndex)    = maxCurrent.at(gaussianIndex).data();
    sigmaCurrentXPtrs.at(gaussianIndex) = sigmaCurrentX.at(gaussianIndex).data();
    sigmaCurrentTPtrs.at(gaussianIndex) = sigmaCurrentT.at(gaussianIndex).data();
  }
  if (infoGlobal->dimensions == 2) {
    calculateDistance = [this](PopInt neuronPop, NeuronInt neuron, int gaussianIndex) {
      return pow(neurons->GetXPosition(neuronPop, neuron) - gaussianPositionX.at(gaussianIndex) * infoGlobal->xAxisLength, 2) +
             pow(neurons->GetYPosition(neuronPop, neuron) - gaussianPositionY.at(gaussianIndex) * infoGlobal->yAxisLength, 2);
    };
  } else {
    calculateDistance = [this](PopInt neuronPop, NeuronInt neuron, int gaussianIndex) {
      return pow(neurons->GetXPosition(neuronPop, neuron) - gaussianPositionX.at(gaussianIndex) * infoGlobal->xAxisLength, 2);
    };
  }
  // Precalc the muMax factor
  for (int gaussianIndex : std::ranges::views::iota(0, noGaussians)) {
    double sigmaX = sigmaCurrentXPtrs.at(gaussianIndex)->parameterValues.at(0);
    for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
      double scalingConstantDt = infoGlobal->dtTimestep * pow(GetScaling(neuronPop), -(infoGlobal->networkScaling_synStrength));
      for (NeuronInt neuron : std::ranges::views::iota(0, neurons->GetNeuronsPop(neuronPop))) {
        // double distance= pow(neurons->GetXPosition(neuronPop, neuron) - gaussianPositionX.at(gaussianIndex)*infoGlobal->xAxisLength, 2) +
        // pow(neurons->GetYPosition(neuronPop, neuron) - gaussianPositionY.at(gaussianIndex)*infoGlobal->yAxisLength, 2);
        cachedMuNeuronFactors.at(gaussianIndex).at(neuronPop).at(neuron) =
            scalingConstantDt * exp(-calculateDistance(neuronPop, neuron, gaussianIndex) / (2 * pow(sigmaX, 2)));
      }
    }
  }
  std::cout << "Done.\n" << std::endl;
}

void SpatialGaussianStimulus::RecalculateFactors(int gaussianIndex, double sigmaX) {
  // std::function<double(PopInt,NeuronInt,int)> distanceCalc;
  // if (infoGlobal->dimensions==2){
  // 	distanceCalc = [this](PopInt neuronPop, NeuronInt neuron, int gaussianIndex){return pow(neurons->GetXPosition(neuronPop, neuron) -
  // gaussianPositionX.at(gaussianIndex)*infoGlobal->xAxisLength, 2) + pow(neurons->GetYPosition(neuronPop, neuron) -
  // gaussianPositionY.at(gaussianIndex)*infoGlobal->yAxisLength, 2);}; } else { 	distanceCalc = [this](PopInt neuronPop, NeuronInt neuron, int
  // gaussianIndex){return pow(neurons->GetXPosition(neuronPop, neuron) - gaussianPositionX.at(gaussianIndex)*infoGlobal->xAxisLength, 2) +
  // pow(neurons->GetYPosition(neuronPop, neuron) - gaussianPositionY.at(gaussianIndex)*infoGlobal->yAxisLength, 2);};
  // }
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    double scalingConstantDt = infoGlobal->dtTimestep * pow(GetScaling(neuronPop), -(infoGlobal->networkScaling_synStrength));
    for (NeuronInt neuron : std::ranges::views::iota(0, neurons->GetNeuronsPop(neuronPop))) {
      // double distance{};
      // if (infoGlobal->dimensions==2){
      // 	distance= pow(neurons->GetXPosition(neuronPop, neuron) - gaussianPositionX.at(gaussianIndex)*infoGlobal->xAxisLength, 2) +
      // pow(neurons->GetYPosition(neuronPop, neuron) - gaussianPositionY.at(gaussianIndex)*infoGlobal->yAxisLength, 2); } else
      // if(infoGlobal->dimensions==1) { 	distance= pow(neurons->GetXPosition(neuronPop, neuron) -
      // gaussianPositionX.at(gaussianIndex)*infoGlobal->xAxisLength, 2);
      // }
      cachedMuNeuronFactors.at(gaussianIndex).at(neuronPop).at(neuron) =
          scalingConstantDt * exp(-calculateDistance(neuronPop, neuron, gaussianIndex) / (2 * pow(sigmaX, 2)));
    }
  }
}

void SpatialGaussianStimulus::SetSignalMatrix() {

  double muMax, sigmaT;
  // double dtTimestep      = infoGlobal->dtTimestep;
  double dtSqrt          = infoGlobal->dtSqrt;
  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  long   timeStep        = infoGlobal->timeStep;
  double mu;
  // double distGaussCenter;//???????

  if ((backgroundNoisePtr != &backgroundNoise.back()) && (backgroundNoisePtr->endTimeStep <= timeStep)) {
    backgroundNoisePtr++;
  }

  for (int gaussianIndex : std::ranges::views::iota(0, noGaussians)) {

    if ((maxCurrentPtrs.at(gaussianIndex) != &maxCurrent.at(gaussianIndex).back()) && (maxCurrentPtrs.at(gaussianIndex)->endTimeStep <= timeStep)) {
      maxCurrentPtrs.at(gaussianIndex)++;
    }
    if ((sigmaCurrentXPtrs.at(gaussianIndex) != &sigmaCurrentX.at(gaussianIndex).back()) &&
        (sigmaCurrentXPtrs.at(gaussianIndex)->endTimeStep <= timeStep)) {
      sigmaCurrentXPtrs.at(gaussianIndex)++;
      RecalculateFactors(gaussianIndex, sigmaCurrentXPtrs.at(gaussianIndex)->parameterValues.at(0));
    }
    if ((sigmaCurrentTPtrs.at(gaussianIndex) != &sigmaCurrentT.at(gaussianIndex).back()) &&
        (sigmaCurrentTPtrs.at(gaussianIndex)->endTimeStep <= timeStep)) {
      sigmaCurrentTPtrs.at(gaussianIndex)++;
    }

    for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
      muMax  = maxCurrentPtrs.at(gaussianIndex)->parameterValues.at(neuronPop);
      sigmaT = sigmaCurrentTPtrs.at(gaussianIndex)->parameterValues.at(neuronPop);

      for (NeuronInt neuron : std::ranges::views::iota(0, neurons->GetNeuronsPop(neuronPop))) {
        mu = muMax * cachedMuNeuronFactors.at(gaussianIndex).at(neuronPop).at(neuron);
        // mu = 0;//For repeating boundaries
        // for (int x_rep = -1;x_rep <= 1;x_rep++) {
        //	for (int y_rep = -1;y_rep <= 1;y_rep++) {
        //		d = pow(neurons->GetX_Pos(neuronPop, i) - infoGlobal->xAxisLength *(GPos_X[Gauss_index] + x_rep), 2) +
        //pow(neurons->GetY_Pos(neuronPop, i) - infoGlobal->yAxisLength *(GPos_Y[Gauss_index] + y_rep), 2); 		if (infoGlobal->Dimensions == 2 || y_rep
        //== 0) 			mu = mu + mu_max * dtTimestep*pow(s, infoGlobal->networkScaling_synStrength)*exp(-d / (2 * pow(sigmax, 2)));
        //	}
        // }
        if (gaussianIndex == 0) {
          signalMatrix.at(neuronPop).at(neuron) = mu * (1 + sigmaT * standardDistribution(generator)) +
                                                  backgroundNoisePtr->parameterValues.at(neuronPop) * dtSqrt * standardDistribution(generator);
        } else {
          signalMatrix.at(neuronPop).at(neuron) += mu * (1 + sigmaT * standardDistribution(generator));
        }
      }
    }
  }
  // Debugging of test9
  //  double accumulator{};
  //  for (std::vector<double>& signalVector : signalMatrix){
  //  	accumulator=std::accumulate(signalVector.begin(), signalVector.end(), accumulator, [](double accumulator, double value){return
  //  accumulator+value;});
  //  }
  //  std::cout<<"signal sum"<<accumulator<<std::endl;
}

void SpatialGaussianStimulus::Update(std::vector<std::vector<double>> &synaptic_dV) {
  if (!wellDefined) {
    return;
  }
  SetSignalMatrix();
  Stimulus::Update(synaptic_dV);
}

double SpatialGaussianStimulus::GetScaling(PopInt neuronPop) const {
  if (infoGlobal->networkScaling_mode == 1) {
    return (neurons->GetTotalNeurons());
  } else if (infoGlobal->networkScaling_mode == 0) {
    return 1.0;
  } else {
    return -1.0;
  }
}
