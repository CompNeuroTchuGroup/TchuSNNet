
//  Created by Pierre Ekelmans on 20/03/2018.
//  Copyright © 2018 Pierre Ekelmans. All rights reserved.
//

#ifndef SpatialGaussianStimulus_hpp
#define SpatialGaussianStimulus_hpp

#include <stdio.h>

#include <iostream>
#include <vector>
#include <random>
#include <limits.h>
#include "Stimulus.hpp"

class SpatialGaussianStimulus : public Stimulus {
protected:
	int noGaussians{1};
	std::normal_distribution<double> standardDistribution{0.0, 1.0};
    std::vector<std::vector<StepStruct>> maxCurrent;
	std::vector<std::vector<StepStruct>> sigmaCurrentT;
	std::vector<std::vector<StepStruct>> sigmaCurrentX;
	std::vector<StepStruct> backgroundNoise;

	std::vector<StepStruct*> maxCurrentPtrs;
	std::vector<StepStruct*> sigmaCurrentTPtrs;
	std::vector<StepStruct*> sigmaCurrentXPtrs;
	StepStruct* backgroundNoisePtr;

	std::vector<double> gaussianPositionX;
	std::vector<double> gaussianPositionY;

	std::vector<std::vector<std::vector<double>>> cachedMuNeuronFactors;

	std::function<double(PopInt,NeuronInt,int)> calculateDistance;
    
    void SetSignalMatrix() override;
	void PostLoadParameters() override;
	void RecalculateFactors(int gaussianIndex, double sigmaX);
public:

	SpatialGaussianStimulus(std::shared_ptr<NeuronPopSample>  neur,std::vector<FileEntry>& parameters,GlobalSimInfo*  infoGlobal);
    ~SpatialGaussianStimulus() override = default;

    std::string GetType() const override{return IDstringSpatialGaussianStimulus;}
    void        Update(std::vector<std::vector<double>>& synaptic_dV) override;

	double GetScaling(PopInt neuronPop) const override;

    void    SaveParameters(std::ofstream& wParameterStream) const override;
    void    LoadParameters(const std::vector<FileEntry>& parameters) override;

};

#endif /* SpatialGaussianStimulus_hpp */
