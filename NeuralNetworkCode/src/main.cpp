/** @brief Simulation of spiking neural network of deterministic excitatory
 * and inhibitory LIF neurons.
 * The program produces data files from your simulation which contain all
 * your specified variables and the simulated population averaged firing rates
 * and the numerical degree of linear summation.
 *
 * @file main.cpp
 *
 * (c) Max-Planck Institute for Brain Research, Frankfurt am Main
 *
 */

#include <iostream>
#include "./DatafileParser.hpp"
#include "./GlobalFunctions.hpp"
#include "NeuralNetwork.hpp"
// All files have been refactored by Antoni Bertolin. If something goes wrong, contact bertolin@uni-bonn.de.

int main(int argc, char *argv[]) {
  // std::string argvString(argv[0]);
  std::string base;
  std::string inputFileAddress;
  bool        Windows         = false;
  std::string pathToInputFile = "";

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
  Windows = true;
#endif

  if (argc >= 2) {
    base = argv[1];
    if (argc >= 3) {
      inputFileAddress = argv[2];

      pathToInputFile = getPathToInputFile(inputFileAddress, Windows);
    } else if (Windows) {
      inputFileAddress = base + "\\Parameters.txt";

      pathToInputFile = base + "\\";
    } else {
      inputFileAddress = base + "/Parameters.txt";

      pathToInputFile = base + "/";
    }
  } else {
    base             = "";
    inputFileAddress = "Parameters.txt";

    pathToInputFile = base;
  }
  std::cout << "Base folder: " + base + "\n";
  std::cout << "Input file : " + inputFileAddress + "\n";

  //************************************************
  // Read Parameter file
  //************************************************
  struct stat buffer;
  if (stat(inputFileAddress.c_str(), &buffer) != 0) {
    std::cout << "*************************\n";
    std::cout << inputFileAddress << " Input file does not exist\n";
    std::cout << "*************************\n";
    return EXIT_FAILURE;
  }

  // array of FileEntry structs
  std::string line;

  // std::vector<std::string>        allParamStringsVector;//,parameterValues;
  // std::string                     readStringLine;//,parameterName,iterateID;
  std::vector<FileEntry>         parameterEntries;
  std::vector<IterableFileEntry> iterate1Entries;
  std::vector<IterableFileEntry> iterate2Entries;

  std::ifstream parameterFileStream(inputFileAddress, std::ios::in);
  while (std::getline(parameterFileStream, line, '\n')) {
    // while (parameterFileStream.getline(line,256)){
    if (line[0] == '#') {
      continue;
    }
    // allParamStringsVector.push_back(line);
    // readStringLine = line;
    // prefix = readStringLine.substr(0, 9);

    // if the entry in the parameter file is prefixed with iterate_1 or iterate_2, push them on iterateX_entries vector
    if (line.find("iterate_1") != std::string::npos) {  //(prefix.compare("iterate_1") == 0) {
      iterate1Entries.push_back(SplitStringToIterableEntry(line));
    } else if (line.find("iterate_2") != std::string::npos) {
      iterate2Entries.push_back(SplitStringToIterableEntry(line));
    } else {
      parameterEntries.push_back(SplitStringToEntry(line));
    }
  }
  parameterFileStream.close();

  // Get path to input parameters file and save it in the parameterEntries
  parameterEntries.push_back(SplitStringToEntry("pathToInputFile  " + pathToInputFile));
  std::for_each(parameterEntries.begin(), parameterEntries.end(), [](FileEntry &entry) { entry.RemoveCommentsInValues(); });
  // parameterEntries.push_back(SplitStringToEntry("nonIterateTitle  " + pathToInputFile));
  // Check for consistency Iterate 1: do all entries have the same lenght?
  if (iterate1Entries.empty()) {
    iterate1Entries.push_back(IterableFileEntry("iterate_1", "placeholder", { "" }));
  }
  // else {
  //     CheckConsistencyOfIterationParameters(iterate1Entries);
  // }

  // Check for consistency Iterate 2: do all entries have the same lenght?
  if (iterate2Entries.empty()) {
    iterate2Entries.push_back(IterableFileEntry("iterate_2", "placeholder", { "" }));
  }
  // else {
  //     CheckConsistencyOfIterationParameters(iterate2Entries);
  // }
  // signed int iterate1_len = static_cast<signed int>();
  // signed int iterate2_len = static_cast<signed int>();
  //******************************************
  //******************************************
  // #pragma omp parallel for collapse(2)
  for (signed int iterate1Index : std::ranges::views::iota(0, MinIterateParameterSize(iterate1Entries))) {
    for (signed int iterate2Index : std::ranges::views::iota(0, MinIterateParameterSize(iterate2Entries))) {
      std::vector<FileEntry> parameterEntriesCopy = parameterEntries;  // Why is this copy created?
      std::cout << "******************************************" << std::endl;
      std::cout << "iterate1 = " << iterate1Index + 1 << " , iterate2 = " << iterate2Index + 1 << std::endl;
      for (FileEntry &parEntry : parameterEntriesCopy) {
        // Set parameters for Iterate 1
        for (IterableFileEntry &parameterEntry1 : iterate1Entries) {
          if (parEntry.parameterName.compare(parameterEntry1.parameterName) == 0) {
            // This loop will allocate the parameters as long as the iterate parameter is consistent with the actual parameter (this
            // requires the parameter be introduced in the params file properly).
            size_t indexMultiplier { IsIterateParamConsistent(parEntry, parameterEntry1) };
            for (size_t index : std::ranges::views::iota(0ull, parEntry.parameterValues.size())) {
              parEntry.parameterValues.at(index) = parameterEntry1.parameterValues.at(indexMultiplier * iterate1Index + index);
            }
            std::cout << " " << parameterEntry1.parameterName << " = "
                      << (parameterEntry1.parameterValues.at(iterate1Index * indexMultiplier)) << std::endl;
            break;
          }
        }
        // Set parameters for Iterate 2
        for (IterableFileEntry &parameterEntry2 : iterate2Entries) {
          if (parEntry.parameterName.compare(parameterEntry2.parameterName) == 0) {
            // This loop will allocate the parameters as long as the iterate parameter is consistent with the actual parameter (this
            // requires the parameter be introduced in the params file properly).
            size_t indexMultiplier { IsIterateParamConsistent(parEntry, parameterEntry2) };
            for (size_t index : std::ranges::views::iota(0ull, parEntry.parameterValues.size())) {
              parEntry.parameterValues.at(index) = parameterEntry2.parameterValues.at(indexMultiplier * iterate2Index + index);
            }
            // parEntry.parameterValues.at(0) = parameterEntry2.parameterValues.at(iterate2Index);
            std::cout << " " << parameterEntry2.parameterName << " = "
                      << (parameterEntry2.parameterValues.at(iterate2Index * indexMultiplier)) << std::endl;
            break;
          }
        }

        if ((parEntry.parameterName.compare("Title") == 0)) {
          parEntry.parameterValues.push_back(parEntry.parameterValues.at(0));
          if (!(iterate1Entries.at(0).parameterName.find("placeholder") != std::string::npos)) {
            parEntry.parameterValues.at(0)
              .append("_it1_" + std::to_string(iterate1Index + 1) + "_")
              .append(
                iterate1Entries.at(0).parameterName,
                iterate1Entries.at(0).parameterName.length() - std::min(static_cast<int>(iterate1Entries.at(0).parameterName.length()), 10),
                std::min(iterate1Entries.at(0).parameterName.length(), static_cast<size_t>(10)))
              .append("_" + iterate1Entries.at(0).parameterValues.at(iterate1Index));
          }
          if (!(iterate2Entries.at(0).parameterName.find("placeholder") != std::string::npos)) {
            parEntry.parameterValues.at(0)
              .append("_it2_" + std::to_string(iterate2Index + 1) + "_")
              .append(
                iterate2Entries.at(0).parameterName,
                iterate2Entries.at(0).parameterName.length() - std::min(static_cast<int>(iterate2Entries.at(0).parameterName.length()), 10),
                std::min(iterate2Entries.at(0).parameterName.length(), static_cast<size_t>(10)))
              .append("_" + iterate2Entries.at(0).parameterValues.at(iterate2Index));
          }
        }
      }
      std::cout << "******************************************" << std::endl;
      try {
        NeuralNetwork neuralNetwork(base, parameterEntriesCopy);
        neuralNetwork.Simulate();
        neuralNetwork.MakeInputCopies(inputFileAddress);
        DatafileParser parser(neuralNetwork.GetRecorderReference());
        // neuralNetwork.~NeuralNetwork(); //This line is doing shennanigans I think
        parser.parse();
      } catch (const char *message) {
        std::cout << message << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  //******************************************
  //******************************************
  std::cout << "\n" << std::endl;
  return EXIT_SUCCESS;
}
