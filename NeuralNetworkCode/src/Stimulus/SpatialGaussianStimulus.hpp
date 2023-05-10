
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

struct signal_step {
    int    end_time = INT_MAX;
	std::vector<double> values;
} ;

class SpatialGaussianStimulus : public Stimulus
{
protected:
    unsigned seed;
	int Ngauss;
    std::valarray<std::vector<signal_step>> maxCurrent;
	std::valarray<std::vector<signal_step>> sigmaCurrent_t;
	std::valarray<std::vector<signal_step>> sigmaCurrent_x;
	std::vector<signal_step> backgroundNoise;
	std::valarray<double> GPos_X;
	std::valarray<double> GPos_Y;
    std::mt19937 generator;
    
    void SetSignalArray();
public:

	SpatialGaussianStimulus(NeuronPopSample * neur,std::vector<std::string> *input,GlobalSimInfo  * info);
    virtual ~SpatialGaussianStimulus()=default;

    std::string GetType(){return str_spatialgaussianStimulus;}
    void        Update(std::vector<std::vector<double>> * synaptic_dV);

	virtual double GetScaling(int pop) {
		if (info->networkScaling_extern == 1)
			return (neurons->GetTotalNeurons());
		else if (info->networkScaling_extern == 0) {
			return 1.0;
		}
        else
            return -1.0;
	}

    void    SaveParameters(std::ofstream * stream);
    void    LoadParameters(std::vector<std::string> *input) override;

};

#endif /* SpatialGaussianStimulus_hpp */
