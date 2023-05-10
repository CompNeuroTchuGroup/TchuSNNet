//
//  WhiteNoiseStimulus.hpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 02/03/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#ifndef WhiteNoiseStimulus_hpp
#define WhiteNoiseStimulus_hpp

#include <stdio.h>

#include <iostream>
#include <vector>
#include <random>
#include <limits.h>
#include "Stimulus.hpp"

struct step {
    int    end_time = INT_MAX;
    std::vector<double> values;
} ;


class WhiteNoiseStimulus : public Stimulus
{
protected:
    unsigned seed;

    std::vector<step> meanCurrent;     // indices of all neurons that have emitted a spike in the previous time step
    std::vector<step> sigmaCurrent;

    std::mt19937 generator;
    
    void SetSignalArray();
public:

    WhiteNoiseStimulus(NeuronPopSample * neur,std::vector<std::string> *input,GlobalSimInfo  * info);
    virtual ~WhiteNoiseStimulus(){}

    std::string GetType(){return str_whitenoiseStimulus;}
    void        Update(std::vector<std::vector<double>> * synaptic_dV);

    virtual double GetScaling(int pop){
        if(info->networkScaling_extern == 0)
            return 1.0;
        else{
            std::cout << "WhiteNoiseStimulus (GetScaling): External Scaling not well defined!";
            std::terminate();
            return 0;}
    };

	
    void    SaveParameters(std::ofstream * stream);
    void    LoadParameters(std::vector<std::string> *input) override;

};

#endif /* WhiteNoiseStimulus_hpp */
