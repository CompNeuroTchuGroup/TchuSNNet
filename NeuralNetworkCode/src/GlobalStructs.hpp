#ifndef _GLOBAL_STRUCTS_HEADER
#define _GLOBAL_STRUCTS_HEADER

#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

using PopInt    = signed int;
using NeuronInt = signed long;
using TStepInt  = long;
using SynInt    = long;
class Synapse;
using SynapsePtr = std::unique_ptr<Synapse>;
class NeuronPop;
using PopPtr  = std::shared_ptr<NeuronPop>;
using cPopPtr = std::shared_ptr<const NeuronPop>;
struct BaseSynapseSpine;
using BaseSpinePtr = BaseSynapseSpine *;
struct BranchedSynapseSpine;
using BranchedSpinePtr = BranchedSynapseSpine *;
struct AlphaSynapseSpine;
using AlphaSpinePtr = AlphaSynapseSpine *;
struct CoopSynapseSpine;
using CoopSpinePtr = std::unique_ptr<CoopSynapseSpine>;
struct MACRbPSynapseSpine;
using MACRbpSpinePtr = MACRbPSynapseSpine *;
struct Branch;
using BranchPtr = Branch *;
struct AlphaBranch;
using AlphaBranchPtr = AlphaBranch *;

struct GlobalSimInfo {
    std::mt19937 globalGenerator;
    std::string  pathToInputFile { "" };
    int          globalSeed { -1 };
    TStepInt     timeStep {};  // long is long int
    // int     waitingIndex{};
    double    dtTimestep {};
    double    dtSqrt {};
    int       density {};
    double    xAxisLength {};  // Always a positive quantity (number of neurons/dimensions)
    double    yAxisLength {};  // Always a positive quantity (number of neurons/dimensions)
    int       dimensions {};
    double    simulationTime {};
    double    networkScaling_synStrength {};
    int       networkScaling_mode {};  // 0: default (=1)
    NeuronInt totalNeurons {};
    bool      isMock { false };
};

struct binData {
    std::vector<double> potential;    // used to compute average potential of neurons [per population]
    std::vector<double> spikerRatio;  // number of spiked neurons [per population]

    std::vector<double>
      externalCurrent;  // average external input current [to population totalNeuronPopulation] //to vector
    std::vector<std::vector<double>>
      synapticCurrents;  // average current [to population P1][from population P2] //to vector
    // std::valarray<double>  totalCurrent_mean;             // mean of the total input current     [to population
    // totalNeuronPopulation]

    std::vector<std::vector<double>> neuronTotalCurrentMean;  // mean of the total input current     [to neuron N]
    std::vector<double>
      totalCurrentSquared_mean;  // mean of the squared total input current to each neuron [averaged over population Ps]
    // std::valarray<int>               spiker_pop; //Member is currently not in use

    /** synaptic statistics per time step for every postsynaptic and every
     * presynaptic population and every synapse specific data column.
     */
    std::vector<std::vector<std::vector<double>>>
      synapticState;  // synaptic data of synapses [to population P1][from population P2][data entry j] //sum
                      // of vectors, valarray is valid (in last level)
    std::vector<std::vector<double>>      heatmap;  // firing rate [of population i][in each pixel]
    std::vector<std::vector<signed long>> noRecordedSynapses;
};
struct FileEntry {
    std::string              parameterName;
    std::vector<std::string> parameterValues;

    FileEntry(std::string parameterName, std::vector<std::string> parameterValues):
        parameterName(std::move(parameterName)), parameterValues(std::move(parameterValues)) { }
    FileEntry() { }

    void RemoveCommentsInValues(char commentCharacter = '#');
    bool parameterNameContains(std::string stringFilter) const {
      return parameterName.find(stringFilter) != std::string::npos;
    }
};

struct IterableFileEntry : FileEntry {
    std::string iterateID;
    IterableFileEntry(std::string iterateID, std::string parameterName, std::vector<std::string> parameterValues):
        FileEntry(parameterName, parameterValues), iterateID(iterateID) { }
};

struct RecorderOpenStreams {
    std::vector<std::ofstream> heatmapStreamVector;
    std::ofstream              averagesFileStream;
    std::ofstream              rasterplotFileStream;
    std::ofstream              potentialFileStream;
    std::ofstream              currentsFileStream;
    std::ofstream              cCurrentsFileStream;
    std::ofstream              histogramFileStream;
    std::ofstream              meanCorrFileStream;
    std::ofstream              pairCorrFileStream;
    std::ofstream              synStatesFileStream;
    std::ofstream              heteroSynFileStream;
    std::ofstream              hSOverallFileStream;
    // std::vector<std::ofstream> neuronOuputFileStreams;
};
struct DendriticSubRegion {  // Still under work
    char             regionID {};
    std::vector<int> branchesInRegion;  // I will have to read this from the morphology LP, every DendriticSubRegion is
                                        // a line, first input is ID, rest is branchIDs. Then in Synapse you put the
                                        // DendriticSubRegion where the synapse goes.
    DendriticSubRegion(char regionID, std::vector<int> branchesInRegion);
};

struct BranchTargeting {  // This is esentially a wrapper for HCS different targeting strategies of postsynaptic
                          // branches
    int              targetBranch {};
    char             DendriticSubRegion { '0' };
    bool             setTargetBranch { false };
    bool             randomTargetBranch { false };
    bool             orderedTargetBranch { false };
    bool             firstSlotTrueLastSlotFalse { true };
    std::vector<int> setOfPositions;
};


struct Instruction {
    // This struct allows the reading organization of instruction
    NeuronInt neuronId {};  // The neuron specific to the instruction (-1 is equivalent to all?)
    long startTimeStep {};  // Will have to convert times to timesteps. Or make a short python programme to do it in a
                            // file by itself
    long   endTimeStep {};
    double frequency {};
    double firingProbability {};
    long   fireEveryNSteps {};  // This variable describes the timesteps between every AP. If 3, it will fire every 3
                                // timesteps, so 2 no and 1 yes.
    bool last { false };
    bool off { false };
    Instruction(FileEntry inputEntry, double dtTimestep);
};

#endif