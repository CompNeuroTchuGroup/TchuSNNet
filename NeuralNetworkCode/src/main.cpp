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
#include "./DatafileParser.hpp"
#include "NeuralNetwork.hpp"
#include "GlobalFunctions.hpp"
#include <iostream>
//#include <filesystem>
//namespace fs = std::__fs::filesystem;

int main(int argc, char* argv[])
{
    std::string argv_str(argv[0]);
    std::string base;
    std::string inputFile;
	bool Windows=false;
    std::string pathTo_inputFile = "";


	#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
		Windows = true;
	#endif

    if(argc >= 2){
        base = argv[1];
        if(argc >= 3){
            inputFile = argv[2];

            pathTo_inputFile = getPathToInputFile(&inputFile, Windows);
        }
        else if (Windows) {
            inputFile = base + "\\Parameters.txt";

            pathTo_inputFile = base + "\\";
        }
        else {
            inputFile = base + "/Parameters.txt";

            pathTo_inputFile = base + "/";
        }
    }
    else{
        base = "";
        inputFile = "Parameters.txt";

        pathTo_inputFile = base;
    }
    std::cout<<"Base folder: "+base+"\n";
    std::cout<<"Input file : "+inputFile+"\n";

    //******************************************
    // Read parameter file
    //******************************************

    //array of FileEntry structs
    char                            line[2048];
    std::vector<FileEntry> parameterEntries;
    std::vector<std::string>        full_strs,values;
    std::string                     str_line,name,iterate_id;
    std::vector<IterableFileEntry> iterate1_entries;
    std::vector<IterableFileEntry> iterate2_entries;

    //************************************************
    // Read Parameter file
    //************************************************
    struct stat buffer;
    if(stat(inputFile.c_str(),&buffer) != 0){
        std::cout << "*************************\n";
        std::cout <<inputFile<<" Input file does not exist\n";
        std::cout << "*************************\n";
        return 0;
    }

    std::ifstream stream(inputFile);
    //std::cout << inputFile << std::endl;

    std::string prefix;
    while (stream.getline(line,256)){
        if(line[0] == '#')
            continue;

        full_strs.push_back(line);
        str_line = line;
        prefix = str_line.substr(0, 9);

        // if the entry in the parameter file is prefixed with iterate_1 or iterate_2, push them on iterateX_entries vector
        if (prefix.compare("iterate_1") == 0) {
            iterate1_entries.push_back(stringToIterableFileEntry(str_line));
        } else if (prefix.compare("iterate_2") == 0) {
            iterate2_entries.push_back(stringToIterableFileEntry(str_line));
        } else {
            parameterEntries.push_back(stringToFileEntry(str_line));
        }
    }
    stream.close();

    //Get path to input parameters file and save it in the parameterEntries
    parameterEntries.push_back(stringToFileEntry("pathToInputFile  " + pathTo_inputFile));

    //Check for consistency Iterate 1: do all entries have the same lenght?
    if(iterate1_entries.empty()){
        iterate1_entries.push_back(IterableFileEntry("iterate_1", "placeholder", {""}));
    } else {
        try {
            checkConsistencyOfIterationParameters(iterate1_entries);
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            abort();
        }
    }

    //Check for consistency Iterate 2: do all entries have the same lenght?
    if(iterate2_entries.empty()){
        iterate2_entries.push_back(IterableFileEntry("iterate_2", "placeholder", {""}));
    } else {
        try {
            checkConsistencyOfIterationParameters(iterate2_entries);
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
            abort();
        }
    }

    unsigned int iterate1_len = static_cast<unsigned int>(iterate1_entries.front().values.size());
    unsigned int iterate2_len = static_cast<unsigned int>(iterate2_entries.front().values.size());

    //******************************************
    //******************************************
    // #pragma omp parallel for collapse(2)
    for(unsigned int i1 = 0; i1 < iterate1_len; i1++) {
        for(unsigned int i2 = 0; i2 < iterate2_len; i2++){

            std::vector<FileEntry> parEntries = parameterEntries;
            std::cout << "******************************************" << std::endl;
            std::cout << "i1 = " << i1+1 << " , i2 = " << i2+1 << std::endl;
            for(auto & parEntry : parEntries){
                //Set parameters for Iterate 1
                for (auto & parEntry_1 : iterate1_entries) {
                    if(parEntry.name.compare(parEntry_1.name) == 0){
                        parEntry.values[0] = parEntry_1.values[i1];
                        std::cout << " " << parEntry_1.name << " = " << (parEntry_1.values[i1]) << std::endl;
                        break;
                    }
                }
                //Set parameters for Iterate 2
                for (auto & parEntry_2 : iterate2_entries) {
                    if(parEntry.name.compare(parEntry_2.name) == 0){
                        parEntry.values[0] = parEntry_2.values[i2];
                        std::cout << " " << parEntry_2.name << " = " << (parEntry_2.values[i2]) << std::endl;
                        break;
                    }
                }

                if((parEntry.name.compare("Title") == 0)){
                    if(iterate1_entries[0].name.compare("placeholder") != 0)
                        parEntry.values[0].append("_i1_"+std::to_string(i1+1)+"_"+iterate1_entries[0].name.substr(iterate1_entries[0].name.length()-6,6)+"_"+iterate1_entries[0].values[i1]);
                    if(iterate2_entries[0].name.compare("placeholder") != 0)
                        parEntry.values[0].append("_i2_"+std::to_string(i2+1)+"_"+iterate2_entries[0].name.substr(iterate2_entries[0].name.length()-6,6)+"_"+iterate2_entries[0].values[i2]);
                }
            }
            std::cout << "******************************************" << std::endl;

            NeuralNetwork neuralNetwork(base,&parEntries);
            neuralNetwork.Simulate();
            neuralNetwork.makeInputCopy(inputFile);
            DatafileParser parser(neuralNetwork.GetRecorder());
            //neuralNetwork.~NeuralNetwork(); //This line is doing shennanigans I think
            parser.parse();
        }
    }


    //******************************************
    //******************************************
    std::cout<< std::endl <<std::endl;
	return 0;
}
