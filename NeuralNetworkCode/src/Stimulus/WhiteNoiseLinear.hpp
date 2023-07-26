//
//  WhiteNoiseLinear.hpp
//  NeuralNetworkCode
//
//  Created by Pierre Ekelmans on 19/02/2018.
//  Copyright © 2018 Pierre Ekelmans. All rights reserved.
//

#ifndef WhiteNoiseLinear_hpp
#define WhiteNoiseLinear_hpp

#include "Stimulus.hpp"
#include <stdio.h>
#include <limits.h>
#include <iostream>
#include <vector>
#include <random>

class WhiteNoiseLinear : public Stimulus {
protected:
	// signed seed{};

	std::vector<StepStruct> meanCurrent;     
	std::vector<StepStruct> sigmaCurrent;

	StepStruct* meanCurrentPtr;
	StepStruct* sigmaCurrentPtr;

	std::vector<double> cachedScalingConstants;
	std::vector<double> meanScalingCache;
	std::vector<double> sigmaScalingCache;

	std::normal_distribution<double> standardDistribution{0.0,1.0};

	// std::mt19937 generator;

	void PostLoadParameters() override;
    void SetSignalMatrix() override;
	void RecalculateMeanScalingCache(StepStruct* meanStep);
	void RecalculateSigmaScalingCache(StepStruct* sigmaStep);
public:

	WhiteNoiseLinear(std::shared_ptr<NeuronPopSample>  neurons, std::vector<FileEntry>& stimulusParameters, GlobalSimInfo*  infoGlobal);
	~WhiteNoiseLinear() override = default;

	std::string GetType() const override { return IDstringWhiteNoiseLinear; }
	void        Update(std::vector<std::vector<double>>& synaptic_dV) override;

	double GetScaling(PopInt neuronPop) const override;

	void    SaveParameters(std::ofstream& wParameterStream) const override;
	void    LoadParameters(const std::vector<FileEntry>& stimulusParameters) override;
};

#endif /* WhiteNoiseLinear_hpp */
