//
//  WhiteNoiseLinear.hpp
//  NeuralNetworkCode
//
//  Created by Pierre Ekelmans on 19/02/2018.
//  Copyright © 2018 Pierre Ekelmans. All rights reserved.
//

#ifndef WhiteNoiseLinear_hpp
#define WhiteNoiseLinear_hpp

#include <stdio.h>

#include <iostream>
#include <vector>
#include <random>
#include <limits.h>
#include "Stimulus.hpp"

struct Step {
	int    start_time = 0;
	int    end_time = INT_MAX;
	std::vector<double> values;
};

class WhiteNoiseLinear : public Stimulus
{
protected:
	unsigned seed;

	std::vector<Step> meanCurrent;     
	std::vector<Step> sigmaCurrent;

	std::mt19937 generator;

	void SetSignalArray();
public:

	WhiteNoiseLinear(NeuronPopSample * neur, std::vector<std::string> *input, GlobalSimInfo  * info);
	virtual ~WhiteNoiseLinear() {}

	std::string GetType() { return str_whitenoiseLinear; }
	void        Update(std::vector<std::vector<double>> * synaptic_dV);

	virtual double GetScaling(int pop) {
		if (info->networkScaling_extern == 0)
			return 1.0;
		else if (info->networkScaling_extern == 1)
			return (neurons->GetTotalNeurons());
		else {
			std::cout << "WhiteNoiseLinear (GetScaling): External Scaling not well defined!";
			std::terminate();
			return 0;
		}
	};

	void    SaveParameters(std::ofstream * stream);
	void    LoadParameters(std::vector<std::string> *input) override;
};

#endif /* WhiteNoiseLinear_hpp */
