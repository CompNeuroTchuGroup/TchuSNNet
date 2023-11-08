//
//  GlobalFunctions.h
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 13/01/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

//#define _NO_PARALLEL_COMPUTING_
#ifdef _NO_PARALLEL_COMPUTING_
#define PAR std::execution::seq
#define PAR_UNSEQ std::execution::unseq
#else
#define PAR std::execution::par
#define PAR_UNSEQ std::execution::par_unseq
#endif

#ifndef GLOBAL_FUNCTIONS_HEADER_
#define GLOBAL_FUNCTIONS_HEADER_
#define _CRT_SECURE_NO_WARNINGS
#include <memory>
#include <random>
#include <ranges>
#include <fstream>
#include <numeric>
#include <deque>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <mutex>
#include <ctime>

//Type alias list

using PopInt = signed int;
using NeuronInt = signed long;
using TStepInt = long;
class Synapse;
using SynapsePtr = std::unique_ptr<Synapse>;
class NeuronPop;
using PopPtr = std::shared_ptr<NeuronPop>;
struct BaseSynapseSpine;
using BaseSpinePtr = BaseSynapseSpine*;
struct BranchedSynapseSpine;
using BranchedSpinePtr = BranchedSynapseSpine*;
struct AlphaSynapseSpine;
using AlphaSpinePtr = AlphaSynapseSpine*;
struct CoopSynapseSpine;
using CoopSpinePtr = std::unique_ptr<CoopSynapseSpine>;
struct CaResSynapseSpine;
using CaRsSpinePtr=CaResSynapseSpine*;
struct Branch;
using BranchPtr = Branch*;
struct AlphaBranch;
using AlphaBranchPtr =AlphaBranch*;

struct GlobalSimInfo {

    std::mt19937 globalGenerator;
    std::string pathToInputFile{""};
    int		globalSeed{-1};
    TStepInt    timeStep{};//long is long int
    // int     waitingIndex{};
    double  dtTimestep{};
    double  dtSqrt{};
    int     density{};
    double  xAxisLength{}; //Always a positive quantity (number of neurons/dimensions)
    double  yAxisLength{}; //Always a positive quantity (number of neurons/dimensions)
    int     dimensions{};
    double  simulationTime{};
    double  networkScaling_synStrength{};
    int     networkScaling_mode{}; // 0: default (=1)
    NeuronInt    totalNeurons{};
    bool isMock{false};
};

struct FileEntry {
    std::string parameterName;
    std::vector<std::string> parameterValues;

    FileEntry(std::string parameterName, std::vector<std::string> parameterValues) : parameterName(std::move(parameterName)), parameterValues(std::move(parameterValues)) {}
    FileEntry(){}

    void RemoveCommentsInValues(char commentCharacter='#');
    bool parameterNameContains(std::string stringFilter) const {return parameterName.find(stringFilter) != std::string::npos;}
};

struct IterableFileEntry : FileEntry {
    std::string iterateID;
    IterableFileEntry(std::string iterateID, std::string parameterName, std::vector<std::string> parameterValues) :  FileEntry(parameterName, parameterValues), iterateID(iterateID){
    }
};



struct RecorderOpenStreams {
    std::vector<std::ofstream> heatmapStreamVector{};
    std::ofstream averagesFileStream;
    std::ofstream rasterplotFileStream;
    std::ofstream potentialFileStream;
    std::ofstream currentsFileStream;
    std::ofstream cCurrentsFileStream;
    std::ofstream histogramFileStream;
    std::ofstream meanCorrFileStream;
    std::ofstream pairCorrFileStream;
    std::ofstream synStatesFileStream;
    std::ofstream heteroSynFileStream;
    std::ofstream hSOverallFileStream;
    //std::vector<std::ofstream> neuronOuputFileStreams;
};
struct DendriticSubRegion{ //Still under work
    const char regionID;
    const std::vector<int> branchesInRegion; //I will have to read this from the morphology LP, every DendriticSubRegion is a line, first input is ID, rest is branchIDs. Then in Synapse you put the DendriticSubRegion where the synapse goes. 
    DendriticSubRegion(char regionID, std::vector<int> branchesInRegion);
};

struct BranchTargeting{ //This is esentially a wrapper for HCS different targeting strategies of postsynaptic branches
    int targetBranch{};
    char DendriticSubRegion{'0'};
    bool setTargetBranch{false};
    bool randomTargetBranch{false};
    bool orderedTargetBranch{false};
    bool firstSlotTrueLastSlotFalse{true};
    std::vector<int> setOfPositions;
};
namespace threadsafe{
    static std::mutex _timeMutex;
    // tm* localtime(const time_t* timer);
    void put_time(time_t timeObj, const char* formatString, std::stringstream& outputString);
}
const std::string IDstringAdjacencyMatrixConnectivity {"AdjacencyMatrixConnectivity"};
const std::string IDstringRandomConnectivity{"RandomConnectivity"};
const std::string IDstringPoissonConnectivity{"PoissonConnectivity"};
const std::string IDstringDistanceConnectivity{"DistanceConnectivity"};
const std::string IDstringHeteroRandomConnectivity{"HeteroRandomConnectivity"};
//const std::string stringlocalConnectivity("LocalConnectivity");


const std::string IDstringCurrentSynapse{"CurrentSynapse"};
const std::string IDstringConductanceSynapse{"ConductanceSynapse"};
const std::string IDstringMongilloSynapse{"MongilloSynapse"};
const std::string IDstringMongilloSynapseContinuous{"MongilloSynapseContinuous"};
const std::string IDstringProbabilisticCurrentSynapse{"ProbabilisticCurrentSynapse"};
const std::string IDstringPRGSynapseContinuous{"PRGSynapseContinuous"};

const std::string IDstringExponentialCurrentSynapse{"ExponentialCurrentSynapse"};
const std::string IDstringPowerLawSynapse{"PowerLawSynapse"};
const std::string IDstringExponentialConductanceSynapse{"ExponentialConductanceSynapse"};
const std::string IDstringExponentialMongilloSynapse{"ExponentialMongilloSynapse"};
const std::string IDstringHeteroSynapse{"HeteroCurrentSynapse"};
const std::string IDstringPlasticityModelSynapse{"PModelSynapse"};

const std::string IDstringExponentialSynapseAddon{"ExponentialSynapseAddon"};

const std::string IDstringLIFNeuron{"LIFNeuron"};
const std::string IDstringQIFNeuron{"QIFNeuron"};
const std::string IDstringEIFNeuron{"EIFNeuron"};
const std::string IDstringPoissonNeuron{"PoissonNeuron"};
const std::string IDstringCorrelatedPoissonNeuron{"CorrelatedPoissonNeuron"};
const std::string IDstringDictatNeuron{"DictatNeuron"};
const std::string IDstringHeteroLIFNeuron{"HeteroLIFNeuron"};
const std::string IDstringHeteroPoissonNeuron{"HeteroPoissonNeuron"};

const std::string IDstringNOPNormalization{"NOPNormalization"};
const std::string IDstringHardNormalization{"HardNormalization"};
const std::string IDstringSoftMaxNormalization{"SoftMaxNormalization"};

//const std::string stringLeanRecorder{"LeanRecorder"};
//const std::string stringAdvancedRecorder{"AdvancedRecorder"};

const std::string IDstringUncorrelatedStimulus{"UncorrelatedStimulus"};
const std::string IDstringWhitenoiseStimulus{"WhiteNoiseStimulus"};
const std::string IDstringWhitenoiseRescaled{"WhiteNoiseRescaled"};
//const std::string stringTims_sin_Stimulus{"TimsSinStimulus"};
const std::string IDstringWhiteNoiseLinear{"WhiteNoiseLinear"};
const std::string IDstringSpatialGaussianStimulus{"SpatialGaussianStimulus"};
const std::string IDstringSpatialPoissonStimulus{"SpatialPoissonStimulus"};

const std::string IDstringMonoDendriteSTDPTazerart{"MonoDendriteSTDPTazerart"};
const std::string IDstringMonoDendriteSTDPTazerartRelative{"MonoDendriteSTDPTazerartRelative"};
const std::string IDstringMonoDendriteSTDPBiWindow{"MonoDendriteSTDPBiWindow"};
const std::string IDstringMonoDendriteSTDPBiExponential{"MonoDendriteSTDPBiExponential"};

const std::string IDstringTraceResourceHSTDP{"AlphaResourceHSTDP"};
const std::string IDstringResourceCalciumDiffusion{"ResourceCalciumDiffusion"};

// void MultiplyVector (std::vector<signed long> &vector, signed long value);
// void MultiplyVector (std::vector<double> &vector, double value);

// void TestWritingFile(std::string filename);

std::vector<FileEntry> FilterStringEntries(const std::vector<FileEntry>& allParameterEntries,std::string filter);
//void FilterStringVector(std::vector<std::string>& fullString,std::string token,std::vector<std::string>& filteredString);
std::vector<std::string> SplitStringToValues(std::string fullString);
FileEntry SplitStringToEntry(std::string fullString);
IterableFileEntry SplitStringToIterableEntry(std::string fullString);
//FileEntry SplitString(std::string& fullString, std::vector<std::string>& parameterValues);

std::string getPathToInputFile(std::string& inputFileAddress, bool Windows);

void SaveDoubleFile(std::ofstream& file,double value,int precision);
// void SaveTupleOfDoublesFile(std::ofstream& file, std::valarray<double>, int precision);
void SaveTupleOfDoublesFile(std::ofstream& file, std::vector<double>, int precision);

bool isDouble(const std::string& readString);

size_t IsIterateParamConsistent(FileEntry entry, IterableFileEntry iterateEntry);
signed int MinIterateParameterSize(std::vector<IterableFileEntry> iterateEntries);
// void CheckConsistencyOfIterationParameters(const std::vector<IterableFileEntry>& iterableEntryVector);

// struct noAllocatableSynapseException : std::exception {
//     char const* what() const noexcept override {
//         return "No synapses available on dendrite for new allocation.";
//     }
// };


//Template functions
//Range()
//V1
//template <typename T>
// auto pyrange(T start, T end, T step) {
//     static_assert(std::is_integral<T>::value, "Integral required.");
//     switch (((start<end) << 1) | (step>0)){
//     case 0:
//         std::swap(start, end);
//         return std::views::iota(0, (end - start + step - 1) / step)
//         | std::views::reverse
//         | std::views::transform([start, step](T value) { return start + value * step; });
//     case 1:
//         std::swap(start, end);
//         return std::views::iota(0, (end - start + step - 1) / step)
//         | std::views::transform([start, step](T value) { return start + value * step; });
//     case 2:
//         continue;
//     case 3:
//         return std::views::iota(0, (end - start + step - 1) / step)
//             | std::views::transform([start, step](T value) { return start + value * step; });
//     case default:
//         throw "Assumptions of ";
//     }
// }
//V2
//Is able to work when compile time evaluatio is not possible
// template <typename T>
// std::vector<T> pyrange(T begin, T end, T stepsize = 1) {
//     std::vector<T> result{};
//     if (begin < end) {
//         for (T item: std::ranges::views::iota(begin, end)
//                  | std::ranges::views::stride(stepsize)) {
//             result.push_back(item);
//         }
//     } else {
//         begin++;
//         end++;
//         stepsize *= -1;
//         for (T item: std::ranges::views::iota(end, begin)         
//                   | std::ranges::views::reverse 
//                   | std::ranges::views::stride(stepsize)) {
//             result.push_back(item);
//         }
//     }
//     return result;
// }
//V3
// template <typename T>
// std::vector<T> pyrange(T begin, T end) {
//     std::vector<T> result{};
//     if (begin < end) {
//         for (T item: std::ranges::views::iota(begin, end)) {
//             result.push_back(item);
//         }
//     } else {
//         begin++;
//         end++;
//         for (T item: std::ranges::views::iota(end, begin)         
//                   | std::ranges::views::reverse ) {
//             result.push_back(item);
//         }
//     }
//     return result;
// }
//V4
//The three dots signal accepting variadic arguments and passing them to iota.
//auto range = [](auto... args) { return std::views::iota(args...); };
//range imitates the behaviour of python's range() with 2 arguments, always with step stride of 1 or -1

#endif /* GlobalFunctions_h */

