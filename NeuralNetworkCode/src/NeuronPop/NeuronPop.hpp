
#ifndef NeuronPop_HPP
#define NeuronPop_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <valarray>
#include <fstream>
#include <string>
#include <memory>
#include "../GlobalFunctions.hpp"
#include "./HeterosynapticNeuronPop/Morphology/SynapseSpines/BaseSynapseSpine.hpp"

typedef std::shared_ptr<BaseSynapseSpine> BaseSpinePtr;

class HeteroCurrentSynapse;
class SynapseSpine;

class NeuronPop
{
protected:

    GlobalSimInfo * info;

    int    identifier;
    unsigned long   noNeurons;
    double refractorySteps;
    double tau_m;
    double v_reset;
    double v_thresh;
    //double r_target;

    int    seed_InitialPotentials,seed_InitialPrevSpike;

    //bool streamingNOutputBool{false};
    bool TaskOutputBool{false};

    std::valarray<double>   potential;         // membrane potential
    std::vector<long>       spiker;            // indices of all neurons that have emitted a spike in the previous time step
    std::valarray<long>     previous_spike_step;   // last spike time for every neuron
	std::valarray<double>   x_pos;   // position on the x axis
	std::valarray<double>   y_pos;   // position on the x axis
    //std::valarray<long>     global_id;  // global neuron id (unique across populations)
    void ClearSpiker();

public:
    // constructor
    NeuronPop(GlobalSimInfo * info,int id);
    virtual ~NeuronPop()=default;

    //*******************
    //Get-Functions
    //*******************
    unsigned long GetNoNeurons(){return noNeurons;}
    //long GetGlobalNeuronId(long i){return global_id[i];}
	double get_Xpos(long i) { return x_pos[i]; }
	double get_Ypos(long i) { return y_pos[i]; }

    // those functions return values and addresses ( pay attention here!!!)
    // TODO/Suggestion: Replace pointer return with a return by constant reference.
    // https://github.com/saiftyfirst/BP_Demos/blob/master/C%2B%2B/constRef_vs_pointer.cpp
    std::vector<long>   * GetSpikers(){return &spiker;}

    //long                get_previous_spike_step(long i) {return previous_spike_step[i];}
    double              GetTimeSinceLastSpike(long i) { return (static_cast<double>(info->time_step) - static_cast<double>(previous_spike_step[i])) * info->dt; }

    double              GetPotential(long i) {return potential[i]; }

    virtual std::string GetType() = 0;
    int  const  GetId(){return this->identifier;}  
    virtual std::valarray<double> GetIndividualSynapticProfile(unsigned long neuronId, unsigned long synapseId);
    virtual std::valarray<double> GetOverallSynapticProfile(unsigned long neuronId);

	//*******************
    //Set-Functions
    //*******************
    void SetNeurons(unsigned long noNeur);
	void SetPosition(long noNeur);
    void SetSeeds(int seed1,int seed2);
    //void SetGlobalNeuronId(long global_neuron_id);

    //Main methods
    virtual void advect(std::vector<double> * synaptic_dV) = 0;
    void    LoadParameters(std::vector<std::string> *input,unsigned long noNeurons);
    virtual void LoadParameters(std::vector<std::string> *input);
    virtual void SaveParameters(std::ofstream * stream);

    //Dyn_casting optimization methods
    virtual bool HasHeterosynapticPlasticity(){return false;}
    virtual bool IsBranchedBool() {return false;}
    virtual unsigned long GetNumberOfSynapses() const {return 0;}// This is purely for virtualization reasons
    //All of the following functions throw to ease the virtualization, but it is a bad coding practice

    virtual void RecordExcitatorySynapticSpike(unsigned long neuronId, unsigned long synapseId);
    virtual BaseSpinePtr AllocateNewSynapse(unsigned long neuronId, HeteroCurrentSynapse& syn);
    virtual std::string GetIndividualSynapticProfileHeaderInfo() const {return "{Something went wrong}";};
    virtual std::string GetOverallSynapticProfileHeaderInfo() const {return "{Something went wrong}";};
    virtual void PostConnectSetUp(){throw;};//This function cannot be called without throwing first
    //Task functions
    //bool const streamingOutput(){return streamingNOutputBool;}//functions
    bool const TaskOutput(){return TaskOutputBool;}
    std::vector<int> TaskOutputVector(); //This function is not complete, will be implemented together with the task framework (supposed to return a vector of 0s and 1s equivalent to spikers)
};


#endif // NeuronPop_HPP
