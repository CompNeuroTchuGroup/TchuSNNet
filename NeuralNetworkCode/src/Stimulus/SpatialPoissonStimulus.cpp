#include "SpatialPoissonStimulus.hpp"

/* When the stimulus has to be changed, the stimulus value at the end
 * is deleted such that at the end of next_stimulus_step there is the current
 * stimulus value.
 * The last value of next_stimulus_time_step is also deleted such that
 * at the end there is the timestep for the next stimulus change.
 */

SpatialPoissonStimulus::SpatialPoissonStimulus(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &stimulusParameters,
                                               GlobalSimInfo *infoGlobal)
    : Stimulus(neurons, infoGlobal) {

  PopInt totalNeuronPops{neurons->GetTotalPopulations()};
  J_External.resize(totalNeuronPops, 0.0);
  lengthScale.resize(totalNeuronPops, 0.0);
  peakConnectProb.resize(totalNeuronPops, 0.0);
  externalCurrents.resize(totalNeuronPops, 0.0);

  LoadParameters(stimulusParameters);

  SetPositions();
  GenerateConnectivity();
  if (inputTau == 0 || inputTau < infoGlobal->dtTimestep) {
    attenuation = 0;
  } else {
    attenuation = 1 - infoGlobal->dtTimestep / inputTau;
  }
}

void SpatialPoissonStimulus::LoadParameters(const std::vector<FileEntry> &stimulusParameters) {

  Stimulus::LoadParameters(stimulusParameters);
  PopInt totalNeuronPops = neurons->GetTotalPopulations();

  for (auto &[parameterName, parameterValues] : stimulusParameters) {

    if ((parameterName.find("virtualExternalNeurons") != std::string::npos) || (parameterName.find("noExternalNeurons") != std::string::npos)) {
      this->noExternalNeurons = std::stoi(parameterValues.at(0));
      sourceFiringPick        = std::uniform_int_distribution<NeuronInt>(0, noExternalNeurons - 1);
    } else if (parameterName.find("stimulus_step") != std::string::npos) {
      AddStimulusStep(std::stod(parameterValues.at(0)), std::stod(parameterValues.at(1)));
    } else if (parameterName.find("J_X") != std::string::npos) {
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        try {
          this->J_External.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
        } catch (...) {
          throw "J_X parameters do not account for all populations";
        }
      }
    } else if (parameterName.find("PeakProba") != std::string::npos) {
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        try {
          this->peakConnectProb.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
        } catch (...) {
          throw "PeakProba parameters do not account for all populations";
        }
      }
    } else if (parameterName.find("lengthscale") != std::string::npos) {
      for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        try {
          this->lengthScale.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
        } catch (...) {
          throw "lengthscale parameters do not account for all populations";
        }
      }
    } else if (parameterName.find("tau_syn") != std::string::npos) {
      this->inputTau = std::stod(parameterValues.at(0));
    } else if (parameterName.find("ExactConnections") != std::string::npos) {
      this->isConnectionExact = std::stoi(parameterValues.at(0));
    }
  }
  if (infoGlobal->isMock) {
    return;
  }
  // With this we can assume the stimulusStep vector will be sorted according to time (reverse sorted before, now because of different finding algo,
  // regular sorted)
  PostLoadParameters();
}

void SpatialPoissonStimulus::SaveParameters(std::ofstream &wParameterStream) const {

  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  Stimulus::SaveParameters(wParameterStream);

  wParameterStream << "stimulus_noExternalNeurons           " << std::to_string(noExternalNeurons) << "\n";
  // if(infoGlobal->globalSeed == -1){
  //     wParameterStream <<  "stimulus_seed                        " << std::to_string(seed)  << "\n";
  // }
  wParameterStream << "stimulus_ExtConnect_lengthscale      ";
  for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
    wParameterStream << std::to_string(lengthScale.at(neuronPop)) << "\t";
  }
  wParameterStream << "#mm\n";
  wParameterStream << "stimulus_ExtConnect_PeakProbability  ";
  for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
    wParameterStream << std::to_string(peakConnectProb.at(neuronPop)) << "\t";
  }
  wParameterStream << "\n";

  wParameterStream << "stimulus_J_X                         ";
  for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
    wParameterStream << std::to_string(J_External.at(neuronPop)) << "\t";
  }
  wParameterStream << "#dmV/Spike\n";

  for (int step : std::ranges::views::iota(0, GetStimulusNoSteps())) {
    wParameterStream << "stimulus_step                        " << std::to_string(GetStimulusStepEndTime(step) * infoGlobal->dtTimestep) << "\t"
                     << std::to_string(GetStimulusStep(step)) << "\t#[t (secs.) -- Hz]\n";
  }

  wParameterStream << "stimulus_tau_syn                     " << std::to_string(inputTau) << "\t\t\t\ts\n";
  wParameterStream << "stimulus_ExactConnections            " << std::to_string(isConnectionExact) << "\t\t\t\t\t\t#(0/1)\n";

  wParameterStream
      << "#\t\t" << IDstringSpatialPoissonStimulus
      << ": noExternalNeurons external neurons with (poisson) firing rate defined in stimulus_step are connected with population X with strength J_X using a distance-dependant connectivity.\n";
}
void SpatialPoissonStimulus::AddStimulusStep(double endTime, double stimStep) {
  StepStruct step;
  step.endTimeStep = static_cast<TStepInt>(std::lround(endTime / infoGlobal->dtTimestep));
  step.parameterValues.push_back(stimStep);
  stimulusSteps.push_back(step);
  // if(stimulusSteps.second.empty())
  // {
  //   stimulusSteps.first.push_back(endTimeStep/infoGlobal->dtTimestep);
  //   stimulusSteps.second.push_back(stimulusStep);
  // }
  // else
  // {
  //   if(endTimeStep/infoGlobal->dtTimestep > stimulusSteps.first.back())
  //   {
  //     std::vector<double> temp_time_step(stimulusSteps.first);
  //     std::vector<double> temp_stimulus(stimulusSteps.second);

  //     stimulusSteps.first.clear();
  //     stimulusSteps.second.clear();

  //     bool done_flag = false;
  //     for(signed i = 0; i < temp_time_step.size(); i++)
  //     {
  //       if(!done_flag && (endTimeStep/infoGlobal->dtTimestep > temp_time_step[i]))
  //       {
  //         stimulusSteps.first.push_back(endTimeStep/infoGlobal->dtTimestep);
  //         stimulusSteps.second.push_back(stimulusStep);
  //         done_flag = true;
  //       }
  //       stimulusSteps.first.push_back(temp_time_step[i]);
  //       stimulusSteps.second.push_back(temp_stimulus[i]);
  //     }
  //     if(stimulusSteps.first.size() == temp_time_step.size())
  //     {
  //       stimulusSteps.first.push_back(endTimeStep/infoGlobal->dtTimestep);
  //       stimulusSteps.second.push_back(stimulusStep);
  //     }
  //   }
  //   else
  //   {
  //     stimulusSteps.first.push_back(endTimeStep/infoGlobal->dtTimestep);
  //     stimulusSteps.second.push_back(stimulusStep);
  //   }
  // }
}

void SpatialPoissonStimulus::PostLoadParameters() {
  std::cout << "\nSetting up the stimulus class...";
  if (stimulusSteps.empty()) { // Nothing should happen if there is no structs in the vector
    wellDefined = false;
    std::cout << "Stimulus class was ill-defined\n" << std::endl;
    return;
  }
  std::sort(stimulusSteps.begin(), stimulusSteps.end(), [](StepStruct &step1, StepStruct &step2) { return step1.endTimeStep < step2.endTimeStep; });
  currentStep = stimulusSteps.data();
  for (PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())) {
    externalCurrents.at(neuronPop) =
        ((J_External.at(neuronPop) * pow(GetScaling(neuronPop), (infoGlobal->networkScaling_synStrength))) * infoGlobal->dtTimestep) / inputTau;
  }
  // Previously in setTableEntries
  //  double signal = infoGlobal->dtTimestep*static_cast<double>(noExternalNeurons)*stimulusSteps.front().parameterValues.at(0);
  poissonDistr =
      std::poisson_distribution<int>(infoGlobal->dtTimestep * static_cast<double>(noExternalNeurons) * stimulusSteps.front().parameterValues.at(0));
  // FillPoissonValueTable(signal);
  std::cout << "Done.\n" << std::endl;
}

void SpatialPoissonStimulus::UpdatePoissonTable() {
  // double  signal;
  // double  dtTimestep = infoGlobal->dtTimestep;
  if ((currentStep != &stimulusSteps.back()) && (currentStep->endTimeStep <= infoGlobal->dtTimestep)) {
    currentStep++;
    poissonDistr =
        std::poisson_distribution<int>(infoGlobal->dtTimestep * static_cast<double>(noExternalNeurons) * currentStep->parameterValues.at(0));
  }
  // StepStruct& step = *std::find_if(stimulusSteps.begin(), stimulusSteps.end(), [this](StepStruct step){
  // 	//Confusing unless you think about a while loop with the condition
  // 	return (infoGlobal->timeStep<=step.endTimeStep) || (step.lastStep);
  // });
  // FillPoissonValueTable(infoGlobal->dtTimestep*static_cast<double>(noExternalNeurons)*step.parameterValues.at(0));
  // if(!stimulusSteps.empty())
  // {
  //    while(stimulusSteps.back().endTimeStep <= infoGlobal->timeStep)
  //     {
  //         stimulusSteps.second.pop_back();
  //         stimulusSteps.first.pop_back();
  //         signal = dtTimestep*static_cast<double>(noExternalNeurons)*stimulusSteps.second.back();
  //         FillPoissonValueTable(signal);

  //         if(stimulusSteps.first.empty())
  //             break;
  //     }
  // }
}
void SpatialPoissonStimulus::SetSignalMatrix() {
  PopInt totalNeuronPops = this->neurons->GetTotalPopulations();
  // int NoFiredExtNeurons;//number of ExternalNeurons that fired during the timestep
  // int sourceNeuronExt;//selected neuron from the external population
  // double current;
  // int targetNeuron;

  UpdatePoissonTable();

  std::for_each(signalMatrix.begin(), signalMatrix.end(), [this](std::vector<double> &signalVector) {
    std::transform(PAR_UNSEQ, signalVector.begin(), signalVector.end(), signalVector.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, this->attenuation));
  });
  // for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)) {
  // 	for (signed long cell : std::ranges::views::iota(0,neurons->GetNeuronsPop(neuronPop))) {
  // 		signalMatrix.at(neuronPop).at(cell) *= attenuation;
  // 	}
  // }

  // NoFiredExtNeurons = static_cast<int>(poissonValueTable.at(tablePick(generator)));
  for (int firedNeuronNo : std::ranges::views::iota(0, poissonDistr(generator))) { // Should this value be rounded?
    (void)firedNeuronNo;
    // For every neuron that fired this timestep: pick the neuron and distribute the current
    NeuronInt sourceNeuronExt = sourceFiringPick(generator);
    for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
      // double current = ((J_External.at(neuronPop) *
      // pow(GetScaling(neuronPop),(infoGlobal->networkScaling_synStrength)))*infoGlobal->dtTimestep)/inputTau;
      std::for_each(externalConnectivity.at(sourceNeuronExt).at(neuronPop).begin(), externalConnectivity.at(sourceNeuronExt).at(neuronPop).end(),
                    [this, neuronPop](NeuronInt synapticTarget) {
                      this->signalMatrix.at(neuronPop).at(synapticTarget) += this->externalCurrents.at(neuronPop);
                    });
      // for (NeuronInt synapticTarget : externalConnectivity.at(sourceNeuronExt).at(neuronPop)) {
      // 	signalMatrix.at(neuronPop).at(synapticTarget) += current;
      // }
    }
  }
}
void SpatialPoissonStimulus::Update(std::vector<std::vector<double>> &synaptic_dV) {
  if (!wellDefined) { // Nothing should happen if there is no structs in the vector
    return;
  }
  SetSignalMatrix();
  Stimulus::Update(synaptic_dV);
}

void SpatialPoissonStimulus::GenerateConnectivity() {
  if (isConnectionExact == 1) {
    GenerateExactConnectivity();
    return;
  }
  double sourcePosX, sourcePosY, targetPosX, targetPosY;
  double distance;
  double probability;
  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  if (infoGlobal->dimensions == 0)
    return;
  externalConnectivity.resize(noExternalNeurons);
  for (NeuronInt sourceExtNeuron = 0; sourceExtNeuron < noExternalNeurons; sourceExtNeuron++) {
    externalConnectivity.at(sourceExtNeuron).resize(totalNeuronPops);
    sourcePosX = externalPosX.at(sourceExtNeuron);
    sourcePosY = externalPosY.at(sourceExtNeuron);
    for (PopInt neuronPop = 0; neuronPop < totalNeuronPops; neuronPop++) {
      for (NeuronInt targetNeuron = 0; targetNeuron < neurons->GetNeuronsPop(neuronPop); targetNeuron++) {
        targetPosX  = neurons->GetXPosition(neuronPop, targetNeuron);
        targetPosY  = neurons->GetYPosition(neuronPop, targetNeuron);
        probability = 0;
        for (int iteratorX = -1; iteratorX <= 1; iteratorX++) {
          if (infoGlobal->dimensions == 1) {
            distance = pow(sourcePosX - targetPosX + iteratorX * infoGlobal->xAxisLength, 2);
            probability += peakConnectProb.at(neuronPop) * exp(-distance / (2 * pow(lengthScale.at(neuronPop), 2)));
          } else if (infoGlobal->dimensions == 2) {
            for (int iteratorY = -1; iteratorY <= 1; iteratorY++) {
              distance = pow(sourcePosX - targetPosX + iteratorX * infoGlobal->xAxisLength, 2) +
                         pow(sourcePosY - targetPosY + iteratorY * infoGlobal->yAxisLength, 2);
              probability += peakConnectProb.at(neuronPop) * exp(-distance / (2 * pow(lengthScale.at(neuronPop), 2)));
            }
          }
        }
        if (zeroToOneDistr(generator) < probability) {
          externalConnectivity.at(sourceExtNeuron).at(neuronPop).push_back(targetNeuron);
        }
      }
    }
  }
}

void SpatialPoissonStimulus::GenerateExactConnectivity() {
  double sourcePosX, sourcePosY, targetPosX, targetPosY;
  long   targetX, targetY;
  double randomR;     // Radius
  double randomTheta; // Angle
  double length;
  PopInt totalNeuronPops = neurons->GetTotalPopulations();
  long   neuronsToConnect;
  long   xAxisSlots;
  long   noTargetNeurons;
  if (infoGlobal->dimensions != 2) {
    throw "Exact SpacePoisson is only available for 2D systems";
  }
  externalConnectivity.resize(noExternalNeurons);
  for (NeuronInt ExtNeuron : std::ranges::views::iota(0, noExternalNeurons)) {
    externalConnectivity.at(ExtNeuron).resize(totalNeuronPops);
    sourcePosX = externalPosX.at(ExtNeuron);
    sourcePosY = externalPosY.at(ExtNeuron);
    for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
      noTargetNeurons  = neurons->GetNeuronsPop(neuronPop);
      length           = lengthScale.at(neuronPop);
      neuronsToConnect = static_cast<long>(peakConnectProb.at(neuronPop) * noTargetNeurons);
      xAxisSlots       = static_cast<long>(round(sqrt(noTargetNeurons)));
      for (int connection : std::ranges::views::iota(0, neuronsToConnect)) {
        (void)connection;
        randomR     = zeroToOneDistr(generator);
        randomTheta = zeroToOneDistr(generator);
        targetPosX  = sourcePosX + length * sqrt(-2 * log(randomR)) * cos(2 * (4 * atan(1)) * randomTheta);
        targetPosY  = sourcePosY + length * sqrt(-2 * log(randomR)) * sin(2 * (4 * atan(1)) * randomTheta);
        while (targetPosX < 0) {
          targetPosX = targetPosX + infoGlobal->xAxisLength;
        }
        while (targetPosY < 0) {
          targetPosY = targetPosY + infoGlobal->xAxisLength; // Shouldn't this be yAxis????
        }
        targetX = static_cast<long>(round(targetPosX / (infoGlobal->xAxisLength / xAxisSlots)));
        targetY = static_cast<long>(round(targetPosY / (infoGlobal->xAxisLength / xAxisSlots)));
        targetX %= xAxisSlots;
        // while (targetX > (xAxisSlots - 1)){//Is this not just modulus???
        // 	targetX = targetX - xAxisSlots;
        // }
        targetY %= xAxisSlots;
        // while (targetY > (xAxisSlots - 1)){//Is this not just modulus???
        // 	targetY = targetY - xAxisSlots;
        // }
        long targetNeuron = targetX + xAxisSlots * targetY;
        externalConnectivity.at(ExtNeuron).at(neuronPop).push_back(targetNeuron);
      }
      std::sort(externalConnectivity.at(ExtNeuron).at(neuronPop).begin(), externalConnectivity.at(ExtNeuron).at(neuronPop).end());
    }
  }
}

// see NeuronPop::SetPosition()
void SpatialPoissonStimulus::SetPositions() {

  int rows;
  int columns;
  int mod;
  externalPosX.resize(noExternalNeurons);
  externalPosY.resize(noExternalNeurons);
  if (infoGlobal->dimensions == 2) {
    rows    = static_cast<int>(ceil(sqrt(noExternalNeurons)));
    columns = static_cast<int>(floor(noExternalNeurons / rows));
    mod     = noExternalNeurons - rows * columns;
    for (NeuronInt neuron : std::ranges::views::iota(0, mod * (columns + 1))) {
      externalPosX.at(neuron) = (infoGlobal->xAxisLength / (columns + 1)) * (neuron % (columns + 1));
      externalPosY.at(neuron) = (infoGlobal->yAxisLength / rows) * floor(neuron / (columns + 1));
    }
    for (NeuronInt neuron : std::ranges::views::iota(mod * (columns + 1), noExternalNeurons)) {
      externalPosX.at(neuron) = (infoGlobal->xAxisLength / (columns)) * ((neuron - mod) % columns);
      externalPosY.at(neuron) = (infoGlobal->yAxisLength / rows) * floor((neuron - mod) / (columns));
    }
  } else if (infoGlobal->dimensions == 1) {
    for (NeuronInt neuron : std::ranges::views::iota(0, noExternalNeurons)) {
      externalPosX.at(neuron) = neuron * infoGlobal->xAxisLength / noExternalNeurons;
    }
  }
}

double SpatialPoissonStimulus::GetScaling(PopInt neuronPop) const {
  double avgConnExtPop; // the number of connections a neuron gets from the external population (on average for the pop)
  if (infoGlobal->networkScaling_mode == 0) {
    if (infoGlobal->dimensions == 1) {
      avgConnExtPop = static_cast<double>(2 * (4 * atan(1)) * noExternalNeurons / (infoGlobal->xAxisLength) * peakConnectProb.at(neuronPop) *
                                          lengthScale.at(neuronPop));
    } else {
      avgConnExtPop = static_cast<double>(2 * (4 * atan(1)) * noExternalNeurons / (infoGlobal->xAxisLength * infoGlobal->yAxisLength) *
                                          peakConnectProb.at(neuronPop) * pow(lengthScale.at(neuronPop), 2));
    }
  } else if (infoGlobal->networkScaling_mode == 1) {
    avgConnExtPop = static_cast<double>(neurons->GetTotalNeurons());
  } else {
    throw "ERROR: Not proper networkScaling in SpatialPoissonStimulus";
  }

  return avgConnExtPop;
}

// fills the poisson_value_table that is a tabular version of the poisson distribution
// with the needed mean that is firing_rate*dt*number_of_neurons
// void SpatialPoissonStimulus::FillPoissonValueTable(double mu){
// 	//What is happening here? Maybe change to usage of poisson distribution?
// 	int value = 0;
// 	double probability = exp(-mu);
// 	double cumulativeProb = exp(-mu);

// 	for (int tableEntry : std::ranges::views::iota(0,tableEntries)) { //Is this a cdf?
// 		if (cumulativeProb < static_cast<double>(tableEntry) / static_cast<double>(tableEntries)) {
// 			value++;
// 			probability *= mu / static_cast<double>(value);
// 			cumulativeProb += probability;
// 		}
// 		poissonValueTable.at(tableEntry) = static_cast<double>(value);
// 	}
// }
