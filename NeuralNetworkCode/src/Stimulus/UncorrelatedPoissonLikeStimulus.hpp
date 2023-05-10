//
//  UncorrelatedPoissonLikeStimulus.h
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 23/01/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#ifndef UncorrelatedPoissonLikeStimulus_h
#define UncorrelatedPoissonLikeStimulus_h

#include <iostream>
#include <vector>
#include <random>
#include "Stimulus.hpp"
#define min_(a,b) (((a)<(b)) ? (a):(b))

/* class UncorrelatedPoissonLikeStimulus is a public Stimulus
 * It operates like a virtual external neuron distribution
 * that fires poisson like and that projects to all totalNeurons neurons
 * in the network. Every neuron in the network receives input from
 * connectedNeurons uncorrelated external neurons.
 * The firing rate of this external population is piecewise constant during
 * the specified time intervals.
 * The number of input spikes for all network neurons for the next time step
 * is evaluated at once and stored in the double * signal_array.
 * Note that not a poisson generator from the standard library is used but
 * a custom made that is less precise due to discretization but much faster for
 * higher firing rates.
 */

class UncorrelatedPoissonLikeStimulus : public Stimulus
{
private:

    unsigned long   noExternalNeurons;
    double          *J_X;

    double * poisson_value_table;    // a table of the poisson distribution for the custom made poisson generator:
    int      table_entries;           // length of poisson_value_table
    int      seed;

    std::vector<double> next_stimulus_time_step;
    std::vector<double>        next_stimulus_step;

    std::mt19937 generator;
    std::uniform_int_distribution<int> distribution;

    inline void UpdatePoissonTable();                // fills the signal_array
    inline void fill_poisson_value_table(double mu); // fills the poisson_value_table

    double GetExternalCouplingStrength(int pop){
        double h;
        if (info->networkScaling_extern == 0)
            h = static_cast<double>(noExternalNeurons);
        else if (info->networkScaling_extern == 1)
            h = static_cast<double>(neurons->GetTotalNeurons());
        else{
            throw "ERROR: GetExternalCouplingStrength";
            h = 0;}

        return (J_X[pop] * pow(h,info->networkScaling_synStrength));

    }

public:
    UncorrelatedPoissonLikeStimulus(NeuronPopSample *neur,std::vector<std::string> *input,GlobalSimInfo  * info);
    ~UncorrelatedPoissonLikeStimulus();

    //*******************
    // Get-Functions
    //*******************
    long    GetStimulusStep_Time(int i) {return static_cast<long> (next_stimulus_time_step.at(i)); }
    double  GetStimulusStep(int i)      {return next_stimulus_step.at(i);};
    long    GetStimulusNoSteps() { return static_cast<long>(next_stimulus_step.size()); }
    std::string GetType()               {return str_uncorrelatedStimulus;}
    int     GetTable_entries()          {return table_entries;}

    //*******************
    // Set-Functions
    //*******************
    void SetSeed(int seed){this->seed = seed; generator = std::mt19937(seed);};
    void SetTableEntries();
    void AddStimulusStep(int ts,double sS);

    //*******************************************
    void    Update(std::vector<std::vector<double>> * synaptic_dV);

    void    SaveParameters(std::ofstream * stream);
    void    LoadParameters(std::vector<std::string> *input) override;
    //void    LoadParameters(std::vector<std::string> *input,double synapticScaling);
};


#endif /* UncorrelatedPoissonLikeStimulus_h */
