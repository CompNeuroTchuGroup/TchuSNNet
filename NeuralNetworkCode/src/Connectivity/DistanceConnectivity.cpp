#include "DistanceConnectivity.hpp"
#include "../Synapse/Synapse.hpp" //NEEDED?

DistanceConnectivity::DistanceConnectivity(Synapse *synapse, GlobalSimInfo *infoGlobal) : Connectivity(synapse, infoGlobal) {
}

void DistanceConnectivity::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
  Connectivity::SaveParameters(wParameterStream, idString);
  wParameterStream << idString << "connectivity_ConnectionProba\t\t\t" << std::to_string(this->peakProbability) << "\n";
  wParameterStream << idString << "connectivity_StdProbability\t\t\t" << std::to_string(this->stdProbability) << " #mm\n";
  wParameterStream << idString << "connectivity_ExactConnections\t\t\t" << std::to_string(this->isConnectionExact)
                   << " #(0/1)\tIf 1, each neuron will receive exactly C connections. \n";
  if (infoGlobal->dimensions == 2) {
    wParameterStream << idString << "connectivity_AsymmetryFactor\t\t\t" << std::to_string(this->projectionLengthFactor)
                     << " # sigma_y=sigma/AF. If AF is 1, the projection length is the same in all directions\n";
  }
  // if (infoGlobal->globalSeed == -1) {
  // 	std::cout << "Please set globalseed >= 0 to use DistanceConnectivity";
  // 	throw "Globalseed is -1";
  // }
  // wParameterStream << "# If ExactConnections==1, neurons have a predefined number of connections randomly generatd from its neighbours . The number
  // of connections is ConnectProba*N.\n"; wParameterStream << "# If ExactConnections==0, the connection probability between each pair of neurons is
  // calculated independantly. There is no repeat, and the number of presynaptic connections is random. ConnectProba defines the peak probability of
  // connection (at zero distance))\n"; wParameterStream << "# "<< stringDistanceConnectivity << ": Each neuronal pair is connected with probability
  // depending on their distance. (as used by [Chariker et al, 2016]). \n";
}

void DistanceConnectivity::LoadParameters(const std::vector<FileEntry> &connectivityParameters) {
  Connectivity::LoadParameters(connectivityParameters);

  for (auto &[parameterName, parameterValues] : connectivityParameters) {
    if (parameterName.find("ConnectionProba") != std::string::npos || parameterName.find("ConnectProba") != std::string::npos) {
      this->peakProbability = std::stod(parameterValues.at(0));
    } else if (parameterName.find("StdProbability") != std::string::npos) {
      this->stdProbability = std::stod(parameterValues.at(0));
    } else if (parameterName.find("ExactConnections") != std::string::npos) {
      this->isConnectionExact = std::stoi(parameterValues.at(0));
    } else if (parameterName.find("AsymmetryFactor") != std::string::npos) {
      this->projectionLengthFactor = std::stoi(parameterValues.at(0));
    }
  }
  // if (infoGlobal->globalSeed == -1) {
  // 	std::cout << "Please set globalseed >= 0 to use DistanceConnectivity";
  // 	throw "Globalseed is -1";
  // }
}

void DistanceConnectivity::ConnectNeurons() {
  if ((infoGlobal->dimensions != 2) && (infoGlobal->dimensions != 1)) {
    std::cout << "Dimensions is not valid for DistanceConnectivity !\n";
    throw "Dimensions is not valid for DistanceConnectivity !\n";
  } else if (isConnectionExact == 1) {
    ConnectNeuronsExact();
    return;
  }
  double      totalConnectionProbability; // probability of connecting neurons Source and Target
  double      randomGeneratedFloat;       // randomly generated number to be compared with totalConnectionProbability
  double      distance;                   // Euclidian distance squared
  NeuronInt   noTargetNeurons = synapse->GetNoTargetNeurons();
  NeuronInt   noSourceNeurons = synapse->GetNoSourceNeurons();
  double      xPositionSource, yPositionSource, xPositionTarget, yPositionTarget;
  NeuronInt   outputInterval = noTargetNeurons / 10;
  signed long synapseCounter = 0;

  if (outputInterval == 0) {
    outputInterval = 1;
  }

  std::uniform_real_distribution<double> distribution(0, 1);

  if (stdProbability > 0.29 * infoGlobal->xAxisLength * std::min(1.0, projectionLengthFactor)) {
    // with std=0.29, the periodic boundaries cause an increase in connection probability of 1% at the center of the domain
    std::cout << "\n \nWarning sigma_ConnectionProbability is high compared to the size of the system \n";
  }
  // Iterate through all target neurons

  for (NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
    // Iterate through all sourceNeuron neurons
    xPositionTarget = synapse->GetTargetNeuronPop()->GetXPosition(targetNeuron);
    yPositionTarget = synapse->GetTargetNeuronPop()->GetYPosition(targetNeuron);

    for (NeuronInt sourceNeuron : std::ranges::views::iota(0, noSourceNeurons)) {
      // Connect with given probability
      xPositionSource = synapse->GetSourceNeuronPop()->GetXPosition(sourceNeuron);
      yPositionSource = synapse->GetSourceNeuronPop()->GetYPosition(sourceNeuron);

      randomGeneratedFloat       = distribution(generator);
      totalConnectionProbability = 0;
      for (int xRepeat :
           std::ranges::views::iota(-1, 2)) { // the motif in which neurons are located is repeated 8 other times to prevent boundary effects
        if (infoGlobal->dimensions == 1) {
          distance = pow((xRepeat * infoGlobal->xAxisLength) + xPositionSource - xPositionTarget, 2); // d^2
          totalConnectionProbability += peakProbability * exp(-distance / (2 * pow(stdProbability, 2)));
        } else if (infoGlobal->dimensions == 2) {
          for (int yRepeat : std::ranges::views::iota(-1, 2)) {
            distance = pow((xRepeat * infoGlobal->xAxisLength) + xPositionSource - xPositionTarget, 2) +
                       pow(projectionLengthFactor, 2) * pow((yRepeat * infoGlobal->yAxisLength) + yPositionSource - yPositionTarget, 2);
            totalConnectionProbability = totalConnectionProbability + peakProbability * exp(-distance / (2 * pow(stdProbability, 2)));
          }
        }
      }
      if ((xPositionTarget == xPositionSource) && (yPositionTarget == yPositionSource) && (sourceNeuron == targetNeuron)) {
        totalConnectionProbability = 0;
      }
      if (randomGeneratedFloat <= totalConnectionProbability) {
        synapse->AllocateSynapse(targetNeuron, sourceNeuron); // synapse->synapticTargets.at(sourceNeuron).push_back(target);
        synapseCounter++;
      }
    }
    // if(targetNeuron%outputInterval == 0)
    //     std::cout << 100*(targetNeuron)/noTargetNeurons << "%\n";
  }
  // std::cout << "100%\n";
  std::cout << " Average Number of Presynaptic neurons : " << std::to_string(synapseCounter / noTargetNeurons) << "\n\n";

  if (stdProbability > 0.29 * infoGlobal->xAxisLength * std::min(1.0, projectionLengthFactor)) {
    std::cout << "Expected number of presynaptic neurons :" << std::to_string(GetExpectedConnections());
    std::cout << " | Actual result : " << std::to_string(synapseCounter / noTargetNeurons) << "\n\n";
  }
}

void DistanceConnectivity::ConnectNeuronsExact() {
  // variation of Connect Neuron, in which the number of neurons connected is set
  // Similar to RandomConnectivity vs BinaryConnectivity
  // can only be used if the number of neurons is a square for each population (sqrt(N_pre) is integer and, sqrt(N_post) is integer)

  NeuronInt   noTargetNeurons = synapse->GetNoTargetNeurons();
  NeuronInt   noSourceNeurons = synapse->GetNoSourceNeurons();
  signed long connectionCounter;
  double      xPositionSource, yPositionSource, xPositionTarget, yPositionTarget;
  NeuronInt   outputInterval = noTargetNeurons / 10;
  long        xSource, ySource;
  NeuronInt   sourceNeuron;
  double      randDistance;
  double      randTheta; // angle
  long        nConnect         = lround((noSourceNeurons - 1) * peakProbability);
  long        xDivisionsSource = noSourceNeurons;                            // NxPre
  double      xDistanceSource  = infoGlobal->xAxisLength / xDivisionsSource; // In one dimension this are the number of divisions //DXsource

  if (infoGlobal->dimensions == 2) {
    // lround() returns a long integer
    xDivisionsSource = lround(sqrt(noSourceNeurons)); // in 2D the distance between rows is different
    xDistanceSource  = infoGlobal->xAxisLength / xDivisionsSource;
    if (0.86 * nConnect > (4 * atan(1)) * pow(2 * stdProbability, 2) * noSourceNeurons / (projectionLengthFactor * pow(infoGlobal->xAxisLength, 2))) {
      // 86% of the Gaussian distribution should fall within the ellipse of radii 2sigma, 2sigma/projectionLengthFactor.
      std::cout
          << "WARNING !! \n Connectivity pattern cannot follow a Gaussian pdf: sigma too low, density too low or Connection_probability too high\n";
    }
  } else if (0.86 * nConnect > 3 * stdProbability * noSourceNeurons / infoGlobal->xAxisLength) {
    // 86% of the Gaussian distribution should fall within the boundary [-1.5sigma 1.5 sigma].
    std::cout
        << "WARNING !! \n Connectivity pattern cannot follow a Gaussian pdf: sigma too low, density too low or Connection_probability too high\n";
  }
  if (pow(xDivisionsSource, 2) != noSourceNeurons && infoGlobal->dimensions == 2) {

    std::cout << "N_pre and N_Post must be exact squares in order to put the neurons on a regular grid in 2D";
    throw "N_pre and N_Post must be exact squares in order to put the neurons on a regular grid in 2D";
  }

  if (outputInterval == 0) {
    outputInterval = 1;
  }

  std::uniform_real_distribution<double> distribution(0, 1);

  for (NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)) {
    // Iterate through all sourceNeuron neurons
    xPositionTarget = synapse->GetTargetNeuronPop()->GetXPosition(targetNeuron);
    yPositionTarget = 0;
    if (infoGlobal->dimensions == 2) {
      yPositionTarget = synapse->GetTargetNeuronPop()->GetYPosition(targetNeuron);
    }
    connectionCounter = 0;
    while (connectionCounter < nConnect) {
      randDistance    = distribution(generator);
      randTheta       = distribution(generator);
      xPositionSource = stdProbability * sqrt(-2 * log(randDistance)) * cos(2 * (4 * atan(1)) * randTheta) + xPositionTarget; // Box Muller transform
      while (xPositionSource < 0) {
        xPositionSource += infoGlobal->xAxisLength;
      }
      xSource = lround(xPositionSource / xDistanceSource);
      xSource %= xDivisionsSource;
      if (infoGlobal->dimensions == 2) {
        yPositionSource =
            stdProbability / projectionLengthFactor * sqrt(-2 * log(randDistance)) * sin(2 * (4 * atan(1)) * randTheta) + yPositionTarget;
        while (yPositionSource < 0) {
          yPositionSource = yPositionSource + infoGlobal->xAxisLength;
        }
        ySource = lround(yPositionSource / xDistanceSource);
        ySource %= xDivisionsSource;
        sourceNeuron = xSource + xDivisionsSource * ySource;
      } else { // 1D{
        sourceNeuron = xSource;
      }
      if ((synapse->IsRecurrent()) && (sourceNeuron == targetNeuron)) {
        continue;
      }
      if (!synapse->IsSourceVectorEmpty(sourceNeuron) && synapse->IsTargetLastInVector(targetNeuron, sourceNeuron)) {
        continue;
      }
      synapse->AllocateSynapse(targetNeuron, sourceNeuron);
      connectionCounter++;
    }

    // if (targetNeuron%outputInterval == 0)
    // 	std::cout << 100 * (targetNeuron) / noTargetNeurons << "%\n";
  }
  // std::cout << "100%\n";
}

double DistanceConnectivity::GetExpectedConnections() const {
  if (isConnectionExact == 1) {
    return (synapse->GetNoSourceNeurons() * peakProbability);
  } else {
    if (infoGlobal->dimensions == 2) {
      return ((synapse->GetNoSourceNeurons() / static_cast<double>(infoGlobal->totalNeurons)) * 2 * (4 * atan(1)) * infoGlobal->density *
              peakProbability * pow(stdProbability, 2) /
              projectionLengthFactor); // integral of density*probability(x,y)dxdy over the infinite 2D space
    } else if (infoGlobal->dimensions == 1) {
      return ((synapse->GetNoSourceNeurons() / static_cast<double>(infoGlobal->totalNeurons)) * sqrt(2 * (4 * atan(1))) * infoGlobal->density *
              peakProbability * stdProbability); // integral of density*probability(x,y)dxdy over the infinite 1D space
    } else {
      throw "ERROR DistanceConnectivity::GetNumberAverageSourceNeurons";
    }
  }
}

// Based on the method of Chariker et al, 2016; and the repeating pattern as in Rosenbaum and Doiron, 2014
