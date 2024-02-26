//
//  GlobalFunctions.h
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 13/01/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

// #define _NO_PARALLEL_COMPUTING_
#ifdef _NO_PARALLEL_COMPUTING_
#define PAR       std::execution::seq
#define PAR_UNSEQ std::execution::unseq
#else
#define PAR       std::execution::par
#define PAR_UNSEQ std::execution::par_unseq
#endif

#ifndef GLOBAL_FUNCTIONS_HEADER_
#define GLOBAL_FUNCTIONS_HEADER_
#define _CRT_SECURE_NO_WARNINGS

#include "GlobalStructs.hpp"

#include <algorithm>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>


// Type alias list

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

namespace threadsafe {
  static std::mutex _timeMutex;
  // tm* localtime(const time_t* timer);
  void put_time(time_t timeObj, const char *formatString, std::stringstream &outputString);
}  // namespace threadsafe
const std::string IDstringAdjacencyMatrixConnectivity { "AdjacencyMatrixConnectivity" };
const std::string IDstringRandomConnectivity { "RandomConnectivity" };
const std::string IDstringPoissonConnectivity { "PoissonConnectivity" };
const std::string IDstringDistanceConnectivity { "DistanceConnectivity" };
const std::string IDstringHeteroRandomConnectivity { "HeteroRandomConnectivity" };
// const std::string stringlocalConnectivity("LocalConnectivity");

const std::string IDstringCurrentSynapse { "CurrentSynapse" };
const std::string IDstringConductanceSynapse { "ConductanceSynapse" };
const std::string IDstringMongilloSynapse { "MongilloSynapse" };
const std::string IDstringMongilloSynapseContinuous { "MongilloSynapseContinuous" };
const std::string IDstringProbabilisticCurrentSynapse { "ProbabilisticCurrentSynapse" };
const std::string IDstringPRGSynapseContinuous { "PRGSynapseContinuous" };

const std::string IDstringExponentialCurrentSynapse { "ExponentialCurrentSynapse" };
const std::string IDstringPowerLawSynapse { "PowerLawSynapse" };
const std::string IDstringExponentialConductanceSynapse { "ExponentialConductanceSynapse" };
const std::string IDstringExponentialMongilloSynapse { "ExponentialMongilloSynapse" };
const std::string IDstringHeteroSynapse { "HeteroCurrentSynapse" };
const std::string IDstringPlasticityModelSynapse { "PModelSynapse" };

const std::string IDstringExponentialSynapseAddon { "ExponentialSynapseAddon" };

const std::string IDstringLIFNeuron { "LIFNeuron" };
const std::string IDstringQIFNeuron { "QIFNeuron" };
const std::string IDstringEIFNeuron { "EIFNeuron" };
const std::string IDstringPoissonNeuron { "PoissonNeuron" };
const std::string IDstringCorrelatedPoissonNeuron { "CorrelatedPoissonNeuron" };
const std::string IDstringDictatNeuron { "DictatNeuron" };
const std::string IDstringHeteroLIFNeuron { "HeteroLIFNeuron" };
const std::string IDstringHeteroPoissonNeuron { "HeteroPoissonNeuron" };

const std::string IDstringNOPNormalization { "NOPNormalization" };
const std::string IDstringHardNormalization { "HardNormalization" };
const std::string IDstringSoftMaxNormalization { "SoftMaxNormalization" };

// const std::string stringLeanRecorder{"LeanRecorder"};
// const std::string stringAdvancedRecorder{"AdvancedRecorder"};

const std::string IDstringUncorrelatedStimulus { "UncorrelatedStimulus" };
const std::string IDstringWhitenoiseStimulus { "WhiteNoiseStimulus" };
const std::string IDstringWhitenoiseRescaled { "WhiteNoiseRescaled" };
// const std::string stringTims_sin_Stimulus{"TimsSinStimulus"};
const std::string IDstringWhiteNoiseLinear { "WhiteNoiseLinear" };
const std::string IDstringSpatialGaussianStimulus { "SpatialGaussianStimulus" };
const std::string IDstringSpatialPoissonStimulus { "SpatialPoissonStimulus" };
const std::string IDstringNoStimulus { "NoStimulus" };

const std::string IDstringMonoDendriteSTDPTazerart { "MonoDendriteSTDPTazerart" };
const std::string IDstringMonoDendriteSTDPTazerartRelative { "MonoDendriteSTDPTazerartRelative" };
const std::string IDstringMonoDendriteSTDPBiWindow { "MonoDendriteSTDPBiWindow" };
const std::string IDstringMonoDendriteSTDPBiExponential { "MonoDendriteSTDPBiExponential" };

const std::string IDstringTraceResourceHSTDP { "AlphaResourceHSTDP" };
const std::string IDstringMACRbPModel { "MACRbPModel" };
const std::string IDstringHeteroGraupnerBrunel { "HeteroGraupnerBrunel" };

// void MultiplyVector (std::vector<signed long> &vector, signed long value);
// void MultiplyVector (std::vector<double> &vector, double value);

// void TestWritingFile(std::string filename);

std::vector<FileEntry> FilterStringEntries(const std::vector<FileEntry> &allParameterEntries, std::string filter);
// void FilterStringVector(std::vector<std::string>& fullString,std::string token,std::vector<std::string>&
// filteredString);
std::vector<std::string> SplitStringToValues(std::string fullString);
FileEntry                SplitStringToEntry(std::string fullString);
IterableFileEntry        SplitStringToIterableEntry(std::string fullString);
// FileEntry SplitString(std::string& fullString, std::vector<std::string>& parameterValues);

std::string getPathToInputFile(std::string &inputFileAddress, bool Windows);

void SaveDoubleFile(std::ofstream &file, double value, int precision);
// void SaveTupleOfDoublesFile(std::ofstream& file, std::valarray<double>, int precision);
void SaveTupleOfDoublesFile(std::ofstream &file, std::vector<double>, int precision);

bool isDouble(const std::string &readString);

size_t     IsIterateParamConsistent(FileEntry entry, IterableFileEntry iterateEntry);
signed int MinIterateParameterSize(std::vector<IterableFileEntry> iterateEntries);
// void CheckConsistencyOfIterationParameters(const std::vector<IterableFileEntry>& iterableEntryVector);

// struct noAllocatableSynapseException : std::exception {
//     char const* what() const noexcept override {
//         return "No synapses available on dendrite for new allocation.";
//     }
// };

// Template functions
// Range()
// V1
// template <typename T>
//  auto pyrange(T start, T end, T step) {
//      static_assert(std::is_integral<T>::value, "Integral required.");
//      switch (((start<end) << 1) | (step>0)){
//      case 0:
//          std::swap(start, end);
//          return std::views::iota(0, (end - start + step - 1) / step)
//          | std::views::reverse
//          | std::views::transform([start, step](T value) { return start + value * step; });
//      case 1:
//          std::swap(start, end);
//          return std::views::iota(0, (end - start + step - 1) / step)
//          | std::views::transform([start, step](T value) { return start + value * step; });
//      case 2:
//          continue;
//      case 3:
//          return std::views::iota(0, (end - start + step - 1) / step)
//              | std::views::transform([start, step](T value) { return start + value * step; });
//      case default:
//          throw "Assumptions of ";
//      }
//  }
// V2
// Is able to work when compile time evaluatio is not possible
//  template <typename T>
//  std::vector<T> pyrange(T begin, T end, T stepsize = 1) {
//      std::vector<T> result{};
//      if (begin < end) {
//          for (T item: std::ranges::views::iota(begin, end)
//                   | std::ranges::views::stride(stepsize)) {
//              result.push_back(item);
//          }
//      } else {
//          begin++;
//          end++;
//          stepsize *= -1;
//          for (T item: std::ranges::views::iota(end, begin)
//                    | std::ranges::views::reverse
//                    | std::ranges::views::stride(stepsize)) {
//              result.push_back(item);
//          }
//      }
//      return result;
//  }
// V3
//  template <typename T>
//  std::vector<T> pyrange(T begin, T end) {
//      std::vector<T> result{};
//      if (begin < end) {
//          for (T item: std::ranges::views::iota(begin, end)) {
//              result.push_back(item);
//          }
//      } else {
//          begin++;
//          end++;
//          for (T item: std::ranges::views::iota(end, begin)
//                    | std::ranges::views::reverse ) {
//              result.push_back(item);
//          }
//      }
//      return result;
//  }
// V4
// The three dots signal accepting variadic arguments and passing them to iota.
// auto range = [](auto... args) { return std::views::iota(args...); };
// range imitates the behaviour of python's range() with 2 arguments, always with step stride of 1 or -1

#endif /* GlobalFunctions_h */
