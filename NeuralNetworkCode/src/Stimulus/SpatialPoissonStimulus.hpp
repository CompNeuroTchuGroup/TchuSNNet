//
//  SpatialPoissonStimulus.h
//  NeuralNetworkCode
//
//  Created by Pierre Ekelmans on 06/04/2018.
//

#ifndef SpatialPoissonStimulus_h
#define SpatialPoissonStimulus_h

#include <iostream>
#include <vector>
#include <random>
#include "Stimulus.hpp"
#define min_(a,b) (((a)<(b)) ? (a):(b))

class SpatialPoissonStimulus : public Stimulus
{
private:
    unsigned long   noExternalNeurons;
	std::vector<double>         J_X;
	std::vector<double>			lengthscale;
	std::vector<double>			PeakConnectProba;
	int exact;

    double * poisson_value_table;    // a table of the poisson distribution for the custom made Poisson generator:
	int      table_entries;           // length of poisson_value_table
    int      seed;
	double tau_syn;
	double attenuation;

	std::valarray<double>   Ext_Xpos;   // position on the x axis
	std::valarray<double>   Ext_Ypos;   // position on the x axis
	std::valarray<std::valarray<std::vector<double>>> Ext_Connectivity;  // Postsynaptic neurons [to external presynaptic neuron i][in population Pi][list of its postsynaptic neurons in this population]
    std::vector<double> next_stimulus_time_step;
    std::vector<double>        next_stimulus_step;

    std::mt19937 generator;
    std::uniform_int_distribution<int> TablePick;
	std::uniform_real_distribution<double> RandomConnect;
	std::uniform_int_distribution<int> FireSourcePick;

    inline void UpdatePoissonTable();                // fills the signal_array
    inline void fill_poisson_value_table(double mu); // fills the poisson_value_table
	void Generate_Connectivity();
	void Generate_Connectivity_Exact();
	void SetPositions();
    double GetExternalCouplingStrength(int pop){
        double h;
		if (info->networkScaling_extern == 0) {
			if (info->Dimensions == 1)
				h = static_cast<double>(2 * (4 * atan(1)) * noExternalNeurons / (info->Lx) * PeakConnectProba[pop] * lengthscale[pop]);
			else
				h = static_cast<double>(2 * (4 * atan(1)) * noExternalNeurons/(info->Lx*info->Ly) * PeakConnectProba[pop] * pow(lengthscale[pop], 2));
		}
        else if (info->networkScaling_extern == 1)
            h = static_cast<double>(neurons->GetTotalNeurons());
        else{
            throw "ERROR: GetExternalCouplingStrength";
            h = 0;}

        return (J_X[pop] * pow(h,info->networkScaling_synStrength));

    }

public:
	SpatialPoissonStimulus(NeuronPopSample *neur,std::vector<std::string> *input,GlobalSimInfo  * info);
    ~SpatialPoissonStimulus();

    //*******************
    // Get-Functions
    //*******************
    long    GetStimulusStep_Time(int i) {return static_cast<long>(next_stimulus_time_step.at(i)); }
    double  GetStimulusStep(int i)      {return next_stimulus_step.at(i);};
    long    GetStimulusNoSteps()        {return static_cast<long> (next_stimulus_step.size());}
    std::string GetType()               {return str_spatialpoissonStimulus;}
    int     GetTable_entries()          {return table_entries;}

    //*******************
    // Set-Functions
    //*******************
    void SetSeed(int seed){this->seed = seed; generator = std::mt19937(seed);};
    void SetTableEntries();
    void AddStimulusStep(double ts,double sS);

    //*******************************************
    void    Update(std::vector<std::vector<double>> * synaptic_dV);

    void    SaveParameters(std::ofstream * stream);
    void    LoadParameters(std::vector<std::string> *input) override;
    //void    LoadParameters(std::vector<std::string> *input,double synapticScaling);
};


#endif /* UncorrelatedPoissonLikeStimulus_h */
