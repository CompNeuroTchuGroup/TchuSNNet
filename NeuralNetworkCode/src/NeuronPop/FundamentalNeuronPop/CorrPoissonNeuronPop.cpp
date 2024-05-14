#include "CorrPoissonNeuronPop.hpp"

CorrPoissonNeuronPop::CorrPoissonNeuronPop(GlobalSimInfo *infoGlobal, NeuronInt neuronID) : NeuronPop(infoGlobal, neuronID) {
  // targetRate = 0; seed = 2;
  //  generator = std::mt19937(seed);
  uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
}

void CorrPoissonNeuronPop::Advect(const std::vector<double> &synaptic_dV) {
  // double dtTimestep           = infoGlobal->dtTimestep;
  ClearSpikerVector();
  if (uniformDistribution(generator) <
      corrLambda) { // this is essentially bernoulli trial, we call Poisson because low prob of 1, negligible of >1 and small timesteps
    spikerNeurons.resize(noNeurons);
    std::iota(spikerNeurons.begin(), spikerNeurons.end(), 0);
  } else {
    NeuronInt totalFiringNeurons{uncorrBinomialDistribution(generator)};
    spikerNeurons.resize(totalFiringNeurons);
    std::sample(neuronIds.begin(), neuronIds.end(), spikerNeurons.begin(), totalFiringNeurons, generator);
  }
  this->AdvectPlasticityModel();
}

void CorrPoissonNeuronPop::LoadParameters(const std::vector<FileEntry> &neuronParameters) {
  NeuronPop::LoadParameters(neuronParameters);

  for (auto &[parameterName, parameterValues] : neuronParameters) {
    if (parameterName.find("r_target") != std::string::npos) {
      finalTargetRate = std::stod(parameterValues.at(0));
      totalLambda     = finalTargetRate * this->infoGlobal->dtTimestep;
      // } else if (parameterName.find("seedPoisson") != std::string::npos) {
      //     seed = std::stoi(parameterValues.at(0));
    } else if (parameterName.find("correlation") != std::string::npos) {
      correlationCoefficient = std::stod(parameterValues.at(0));
      if (correlationCoefficient > 1 || correlationCoefficient < 0) {
        throw "Correlation cannot be outside of the [0,1] range";
      }
      // } else if (parameterName.find("seedPoisson") != std::string::npos) {
      //     seed = std::stoi(parameterValues.at(0));
    }
  }
  PreCalcLambdas();
  // if(infoGlobal->globalSeed != -1){
  //     std::uniform_int_distribution<int> distribution(0,INT32_MAX);
  //     seed = distribution(infoGlobal->globalGenerator);
  //     generator = std::mt19937(seed);
  // }
  neuronIds.resize(noNeurons);
  std::iota(neuronIds.begin(), neuronIds.end(), 0);
  uncorrBinomialDistribution = std::binomial_distribution<>(noNeurons, uncorrLambda);
}

void CorrPoissonNeuronPop::PreCalcLambdas() {
  corrLambda   = correlationCoefficient * totalLambda;
  uncorrLambda = totalLambda - corrLambda;
}

void CorrPoissonNeuronPop::SaveParameters(std::ofstream &wParameterStream) const {

  std::string idString = "neurons_" + std::to_string(GetId());

  // std::cout<<"printing neuron id: "<<id<<"\n";

  // NeuronPop::SaveParameters(stream);
  wParameterStream << "#***********************************************\n";
  wParameterStream << idString + "_noNeurons\t\t\t" << noNeurons << "\n";
  wParameterStream << idString + "_type\t\t\t" << GetType() << "\n";
  wParameterStream << idString + "_correlation_coefficient\t" << std::to_string(correlationCoefficient) << "\tRange from 0 to 1"
                   << "\n";
  // if (infoGlobal->globalSeed == -1) {
  // 	*stream << id + "_seedPoisson                 " << std::to_string(seed) << "\n";
  // }
  wParameterStream << idString + "_r_target\t\t\t" << std::to_string(finalTargetRate) << "\n";

  if (userSeeds) {
    // wParameterStream <<  neuronID + "_seedInitialPotentials   " << this->seedInitialPotentials << "\n";
    // wParameterStream <<  neuronID + "_seedInitialPrevSpike    " << this->seedInitialPreviousSpike << "\n";
    wParameterStream << idString + "_seed\t\t\t\t" << this->seed << "\n";
  }
  wParameterStream
      << "#\t\tPoisson neuron: produces Poisson spiking with rate r_target (defined under stimulus). ZERO DOES NOT REMOVE THE FEATURE, YOU MUST REMOVE THE ENTIRE LINE \n";
  wParameterStream
      << "#\t\tIf r_target not set in parameters, the neurons will fire with probability equal to the membrane potential. If Vm > 1mV, p=1 \n";
}
