#ifndef _DICTAT_NEURON_POP_HEADER_
#define _DICTAT_NEURON_POP_HEADER_

#include "../NeuronPop.hpp"

struct Instruction{
    //This struct allows the reading organization of instruction
    int neuronId;//The neuron specific to the instruction (-1 is equivalent to all?)
    long startTimeStep;//Will have to convert times to timesteps. Or make a short python programme to do it in a file by itself
    long endTimeStep;
    double frequency;
    double lambdaPoisson;
    long fireEveryNSteps;//This variable describes the timesteps between every AP. If 3, it will fire every 3 timesteps, so 2 no and 1 yes.
    bool completed{false};
    bool last{false};
    bool off{false};
    Instruction(int neuronId, double startTime, double endTime, double frequency, double dt);
};

//ControlledNeuronPop
//Reader/FileNeuronPop
//DictatNeuronPop WINNER

class DictatNeuronPop: public NeuronPop {
    protected:
    std::vector<std::vector<Instruction>> inputInstructions;//A vector of instructions for each neuron of the population
    std::vector<int> activeInstructions;
    bool instructionBasedBool{false};
    bool poissonLikeFiring{false};

    std::mt19937 randomGenerator;
    std::uniform_real_distribution<double> uniformDistribution = std::uniform_real_distribution<double>(0.0,1.0);

    std::string inputFileAddress;
    std::ifstream inputStreamFile;//A stream to read the spikers per time-step
    bool spikerListFiringBasedBool{false};

    public:
    DictatNeuronPop(GlobalSimInfo* info, int id);
    ~DictatNeuronPop()=default;

    void LoadParameters(std::vector<std::string> *input) override;
    void SaveParameters(std::ofstream * stream) override;
    
    void advect(std::vector<double> * synaptic_dV) override; //Here I will have to set the spikers depending on the instructions (using modular with fireEveryNsteps)
    //Possibility of instead of instructions, giving the spikers per timestep (probably calculate the raw firing first, then the spikers per timestep) and set them in the advect function
    void ReadInstructionsFromFile();
    void GenerateRegularSpikersFromInstructions();
    void GeneratePoissonSpikersFromInstructions();
    void ReadSpikersFromFile();

    std::string GetType() override {return str_DictatNeuron;}
}; 

#endif