//
//  WhiteNoiseStimulus.cpp
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 02/03/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#include "WhiteNoiseStimulus.hpp"


WhiteNoiseStimulus::WhiteNoiseStimulus(NeuronPopSample *neur,std::vector<std::string> *input,GlobalSimInfo  * info):Stimulus(neur,info){

    //******************************
    //***** Default values *********
    //******************************
    seed         = 0;
    //******************************

    LoadParameters(input);
    generator  = std::mt19937(seed);
}


void WhiteNoiseStimulus::LoadParameters(std::vector<std::string> *input){

    //Stimulus::LoadParameters(input);
    std::string              name,token;
    std::vector<std::string> values;
    int                      P = neurons->GetTotalPopulations();
    step s;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);

        if((name.find("seed") != std::string::npos)){
            seed = static_cast<unsigned int>(std::stod(values.at(0)));
        }
        else if((name.find("meanCurrent") != std::string::npos)){
            s.values.clear();
            for(int i = 0;i< P;i++)
                s.values.push_back(std::stod(values.at(i)));
            if(is_double(values.at(P)))
                s.end_time = static_cast<int>(std::round(std::stod(values.at(P)) / info->dt));
            else
                s.end_time = INT_MAX;
			if ((s.end_time < 0) || (s.end_time >static_cast<int>(info->simulationTime / info->dt)))
                s.end_time = static_cast<int>(std::round(info->simulationTime / info->dt));

            meanCurrent.push_back(s);
        }
        else if((name.find("sigmaCurrent") != std::string::npos)){
            s.values.clear();
            for(int i = 0;i< P;i++)
                s.values.push_back(std::stod(values.at(i)));

            if(is_double(values.at(P)))
                s.end_time = static_cast<int>(std::round(std::stod(values.at(P)) / info->dt));
            else
                s.end_time = INT_MAX;

            if((s.end_time < 0) || (s.end_time > std::round(info->simulationTime/info->dt)))
                s.end_time = static_cast<int>(std::round(info->simulationTime / info->dt));

            sigmaCurrent.push_back(s);
        }
    }

    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        seed = distribution(info->globalGenerator);
    }
}


void WhiteNoiseStimulus::SaveParameters(std::ofstream * stream){

    int P        = neurons->GetTotalPopulations();
    Stimulus::SaveParameters(stream);

    if(info->globalSeed == -1){
        *stream <<  "stimulus_seed                        " << std::to_string(seed)  << "\n";
    }

    for(auto const &s : meanCurrent){
        *stream <<  "stimulus_meanCurrent                 ";
        for(int i = 0;i<P;i++)
            *stream << std::to_string(s.values.at(i)) << "\t ";
        *stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
        *stream << " [column 1: input for population 1, column 2: input for pop. 2, ... , last column: time until which input is set. Dimensions: [mV/sec , secs.]\n";
        //*stream << " [mV/sec -- sec]\n";
    }

    for(auto const &s : sigmaCurrent){
        *stream <<  "stimulus_sigmaCurrent                ";
        for(int i = 0;i<P;i++)
            *stream << std::to_string(s.values.at(i)) << "\t ";
        *stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
        *stream << " [column 1: input for population 1, column 2: input for pop. 2, ... , last column: time until which input is set. Dimensions: [mV/sqrt(sec) , secs.]\n";
    }

    if(GetType().compare("WhiteNoiseRescaled")==0)
        *stream <<  "#\t\tRI_{i,ext}/tau_m*dt = meanCurrent_i*dt + sqrt(dt)*sigmaCurrent_i*NormalDistribution(0,1)\n";
}


void WhiteNoiseStimulus::SetSignalArray(){

    double mean,sigma,s;
    double dt      = info->dt;
    int    P       = neurons->GetTotalPopulations();
    double sqrt_dt = sqrt(dt);
    long   t_step  = info->time_step;
    std::normal_distribution<double> distribution(0,1);
    step   *mean_step_current   = &meanCurrent.at(0);
    step   *sigma_step_current  = &sigmaCurrent.at(0);

    while((t_step > mean_step_current->end_time) && (mean_step_current != &meanCurrent.back()))
        mean_step_current++;

    while((t_step > sigma_step_current->end_time) && (sigma_step_current != &sigmaCurrent.back()))
        sigma_step_current++;

    for(int pop = 0;pop<P;pop++){
        mean  = mean_step_current->values.at(pop);
        sigma = sigma_step_current->values.at(pop);
        s = GetScaling(pop);
        for (unsigned long i = 0;i<neurons->GetNeuronsPop(pop); i++)
            signal_array[pop][i] = mean*dt*pow(s,-(info->networkScaling_synStrength)) + sqrt_dt*sigma*distribution(generator);
    }

}

void WhiteNoiseStimulus::Update(std::vector<std::vector<double>> * synaptic_dV)
{
    SetSignalArray();
    Stimulus::Update(synaptic_dV);
}
