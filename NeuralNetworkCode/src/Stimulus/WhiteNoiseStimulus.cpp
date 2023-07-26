//
//  WhiteNoiseStimulus.cpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 02/03/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#include "WhiteNoiseStimulus.hpp"


WhiteNoiseStimulus::WhiteNoiseStimulus(std::shared_ptr<NeuronPopSample> neurons,std::vector<FileEntry>& stimulusParameters,GlobalSimInfo*  infoGlobal): Stimulus(neurons,infoGlobal){

    //******************************
    //***** Default parameterValues *********
    //******************************
    //******************************
    cachedMeans.resize(neurons->GetTotalPopulations());
    cachedSigmas.resize(neurons->GetTotalPopulations());
    LoadParameters(stimulusParameters);
    // generator  = std::mt19937(seed);
}


void WhiteNoiseStimulus::LoadParameters(const std::vector<FileEntry>& stimulusParameters){

    Stimulus::LoadParameters(stimulusParameters);
    PopInt                      totalNeuronPops = neurons->GetTotalPopulations();

    for(auto& [parameterName, parameterValues] : stimulusParameters) {
        // if((parameterName.find("seed") != std::string::npos)){
        //     seed = static_cast<signed int>(std::stod(parameterValues.at(0)));
        if((parameterName.find("meanCurrent") != std::string::npos)){
            StepStruct step;
            for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
                step.parameterValues.push_back(std::stod(parameterValues.at(neuronPop)));
            }
            if(isDouble(parameterValues.at(totalNeuronPops))){
                step.endTimeStep = static_cast<TStepInt>(std::round(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
                if ((step.endTimeStep < 0) || (step.endTimeStep >static_cast<TStepInt>(infoGlobal->simulationTime / infoGlobal->dtTimestep))){
                    step.endTimeStep = static_cast<TStepInt>(std::round(infoGlobal->simulationTime / infoGlobal->dtTimestep));
                }
            } else {
                step.endTimeStep = static_cast<TStepInt>(std::round(infoGlobal->simulationTime / infoGlobal->dtTimestep));;
            }
			
            meanCurrent.push_back(step);
        }
        else if((parameterName.find("sigmaCurrent") != std::string::npos)){
            StepStruct step;
            for(PopInt neuronPop : std::ranges::views::iota (0, totalNeuronPops)){
                step.parameterValues.push_back(std::stod(parameterValues.at(neuronPop)));
            }
            if(isDouble(parameterValues.at(totalNeuronPops))){
                step.endTimeStep = static_cast<int>(std::round(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
                if((step.endTimeStep < 0) || (step.endTimeStep > std::round(infoGlobal->simulationTime/infoGlobal->dtTimestep))){
                    step.endTimeStep = static_cast<int>(std::round(infoGlobal->simulationTime / infoGlobal->dtTimestep));
                }
            } else {
                step.endTimeStep = static_cast<TStepInt>(std::round(infoGlobal->simulationTime / infoGlobal->dtTimestep));
            }
            sigmaCurrent.push_back(step);
        }
    }
	if(infoGlobal->isMock){
		return;
	}
    PostLoadParameters();
}


void WhiteNoiseStimulus::SaveParameters(std::ofstream& wParameterStream) const{

    PopInt totalNeuronPops        = neurons->GetTotalPopulations();
    Stimulus::SaveParameters(wParameterStream);

    // if(infoGlobal->globalSeed == -1){
    //     wParameterStream <<  "stimulus_seed                        " << std::to_string(seed)  << "\n";
    // }

    for(const StepStruct &step : meanCurrent){
        wParameterStream <<  "stimulus_meanCurrent                 ";
		for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
            wParameterStream << std::to_string(step.parameterValues.at(neuronPop)) << "\t ";
        }
        wParameterStream << std::to_string(static_cast<double>(step.endTimeStep)*infoGlobal->dtTimestep) << " \t";
        wParameterStream << " #[column 1: input for population 1, column 2: input for neuronPop. 2, ... , last column: time until which input is set. Dimensions: [mV/sec , secs.]\n";
        //*stream << " [mV/sec -- sec]\n";
    }

    for(const StepStruct &step : sigmaCurrent){
        wParameterStream <<  "stimulus_sigmaCurrent                ";
		for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
            wParameterStream << std::to_string(step.parameterValues.at(neuronPop)) << "\t ";
        }
        wParameterStream << std::to_string(static_cast<double>(step.endTimeStep)*infoGlobal->dtTimestep) << " \t";
        wParameterStream << " #[column 1: input for population 1, column 2: input for neuronPop. 2, ... , last column: time until which input is set. Dimensions: [mV/sqrt(sec) , secs.]\n";
    }

    if(infoGlobal->networkScaling_mode == 1){
        wParameterStream <<  "#\t\tRI_{i,ext}/tauM*dt = meanCurrent_i*dt + sqrt(dt)*sigmaCurrent_i*NormalDistribution(0,1)\n";
        wParameterStream<< "#\t\tmeanCurrent_i is rescaled with N^(-scalingSynapticStrength)";
    }
}

void WhiteNoiseStimulus::RecalculateMeansCache(StepStruct *meanStep) {
    double dtTimestep{infoGlobal->dtTimestep};
	for (PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
		cachedMeans.at(neuronPop) = meanStep->parameterValues.at(neuronPop)*dtTimestep*pow(GetScaling(neuronPop),-(infoGlobal->networkScaling_synStrength));
	}
}

void WhiteNoiseStimulus::RecalculateSigmasCache(StepStruct *sigmaStep) {
    double dtSqrtd{infoGlobal->dtSqrt};
	for (PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
		cachedSigmas.at(neuronPop) = dtSqrtd*sigmaStep->parameterValues.at(neuronPop);
	}
}

void WhiteNoiseStimulus::PostLoadParameters() {
    std::cout<<"\nSetting up the stimulus class...";
    if(meanCurrent.empty() || sigmaCurrent.empty()){
        wellDefined=false;
		std::cout<<"Stimulus class was ill-defined\n"<<std::endl;
        return;
    }

    double dtTimestep{infoGlobal->dtTimestep};
    double dtSqrtd{infoGlobal->dtSqrt};

    currentMeanStep=meanCurrent.data();
    currentSigmaStep=sigmaCurrent.data();

    for (PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
		cachedMeans.at(neuronPop) = currentMeanStep->parameterValues.at(neuronPop)*dtTimestep*pow(GetScaling(neuronPop),-(infoGlobal->networkScaling_synStrength));
	}
    for (PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
		cachedSigmas.at(neuronPop) = dtSqrtd*currentSigmaStep->parameterValues.at(neuronPop);
	}
    std::cout<<"Done.\n"<<std::endl;
}

void WhiteNoiseStimulus::SetSignalMatrix() {

    TStepInt   timeStep = infoGlobal->timeStep;

    // StepStruct   meanStepCurrent   = meanCurrent.at(0);
    // StepStruct   stepCurrentSigma  = sigmaCurrent.at(0);

    //Change this with 
    if((currentMeanStep != &meanCurrent.back()) && (currentMeanStep->endTimeStep<=timeStep)){
		currentMeanStep++;
		RecalculateMeansCache(currentMeanStep);
	}
	// const StepStruct& meanStepCurrent = *std::find_if(meanCurrent.begin(), meanCurrent.end(), [timeStep](StepStruct step){
	// //Same as the while statement
	// 		return (timeStep<=step.endTimeStep) || (step.lastStep);
	// });
    if((currentSigmaStep != &sigmaCurrent.back()) && (currentSigmaStep->endTimeStep<=timeStep)){
		currentSigmaStep++;
		RecalculateSigmasCache(currentSigmaStep);
	}
	// const StepStruct& stepCurrentSigma = *std::find_if(sigmaCurrent.begin(), sigmaCurrent.end(), [timeStep](StepStruct step){
	// //Same as the while statement
	// 		return (timeStep<=step.endTimeStep) || (step.lastStep);
	// });

    for(PopInt neuronPop : std::ranges::views::iota (0, neurons->GetTotalPopulations())){
        for (NeuronInt neuron : std::ranges::views::iota(0,neurons->GetNeuronsPop(neuronPop))){
            signalMatrix.at(neuronPop).at(neuron) = cachedMeans.at(neuronPop) + cachedSigmas.at(neuronPop)*standardDistribution(generator);
        }
    }
}

void WhiteNoiseStimulus::Update(std::vector<std::vector<double>>& synaptic_dV) {
    SetSignalMatrix();
    Stimulus::Update(synaptic_dV);
}

double WhiteNoiseStimulus::GetScaling(PopInt neuronPop) const {
    if(infoGlobal->networkScaling_mode == 0){
        return 1.0;
    } else if(infoGlobal->networkScaling_mode == 1){
        return (neurons->GetTotalNeurons());
    }else{
        throw "WhiteNoiseStimulus (GetScaling): External Scaling not well defined!";
    }
}

