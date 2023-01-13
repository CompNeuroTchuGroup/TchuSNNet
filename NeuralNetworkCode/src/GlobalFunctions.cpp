//
//  GlobalFunctions.c
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 13/01/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <string>


#include "GlobalFunctions.hpp"



void multiply_vector (std::vector<unsigned long> &vector, unsigned long value)
{
    for(std::vector<unsigned long>::iterator it = vector.begin(); it != vector.end(); ++it)
    {
        *it = *it*value;
    }
}


void multiply_vector (std::vector<double> &vector, double value)
{
    for(std::vector<double>::iterator it = vector.begin(); it != vector.end(); ++it)
    {
        *it = *it*value;
    }
}


void TestWritingFile(std::string filename)
{
    // if this test file does not appear in the target directory: stop the
    // simulation and check the directoryPath.
    std::ofstream testStream(filename);
    testStream
    << "Hello World!\n"
    << "This is a serious Panda test!\n"
    << "You passed :-)\n"
    << "Your directory path was correct.\n"
    << std::endl;
    testStream.close();
}


void FilterStringVector(std::vector<std::string> *str_full,std::string token,std::vector<std::string> *str_filtered){
    str_filtered->clear();
    for(std::vector<std::string>::iterator it = (*str_full).begin(); it != (*str_full).end(); ++it) {
        if(it->find(token) != std::string::npos){
            str_filtered->push_back(*it);
        }
    }
}

void FilterStringEntries(std::vector<ParameterFileEntry> *str_full,std::string token,std::vector<ParameterFileEntry> *str_filtered){
    str_filtered->clear();

    for (auto & parEntry : *str_full) {
        if(parEntry.name.find(token) != std::string::npos){
            str_filtered->push_back(parEntry);
        }
    }


}

void SplitString(std::string* full_str, std::vector<std::string>* values) {

    std::replace(full_str->begin(), full_str->end(), '\t', ' '); //Replace tabs by spaces, then search for tabs
    std::replace(full_str->begin(), full_str->end(), '\r', ' ');  //Replace additional end-of-line \r-symbol

    std::istringstream iss(*full_str);
    std::string        token;
    values->clear();

    while (std::getline(iss, token, ' ')) {
        if (token.length() != 0) {
            values->push_back(token);
        }
    }

}

std::string getPathToInputFile(std::string* inputFile, bool Windows) {

    std::istringstream iss(*inputFile);
    std::string        token;
    std::string        pathTo_inpuFile = "";
    char tokenizerChar = ' ';
    if (Windows)
        tokenizerChar = '\\';
    else
        tokenizerChar = '/';

    while (std::getline(iss, token, tokenizerChar)) {
        if (!iss.eof()) {
            if (token.length() != 0) {
                pathTo_inpuFile += token;
                if (Windows)
                    pathTo_inpuFile += "\\";
                else
                    pathTo_inpuFile += "/";
            }
        }
    }

    return pathTo_inpuFile;
}

void SplitString(std::string *full_str,std::string *name,std::vector<std::string> *values){

    std::replace(full_str->begin(),full_str->end(),'\t',' '); //Replace tabs by spaces, then search for tabs
    std::replace(full_str->begin(),full_str->end(),'\r',' ');  //Replace additional end-of-line \r-symbol

    std::istringstream iss(*full_str);
    std::string        token;
    values->clear();

    std::getline(iss,*name,' '); //ignore first entry
    while(std::getline(iss,token,' ')){
        if(token.length() != 0){
            values->push_back(token);
        }
    }

}

void SplitString(std::string *full_str,std::string* iterate_id, std::string *name,std::vector<std::string> *values){

    std::replace(full_str->begin(),full_str->end(),'\t',' '); //Replace tabs by spaces, then search for tabs
    std::replace(full_str->begin(),full_str->end(),'\r',' ');  //Replace additional end-of-line \r-symbol

    std::istringstream iss(*full_str);
    std::string        token;
    values->clear();

    std::getline(iss,*iterate_id,' '); //ignore first two entries
    std::getline(iss,*name,' ');
    while(std::getline(iss,token,' ')){
        if(token.length() != 0){
            values->push_back(token);
        }
    }

}

void SaveDoubleFile(std::ofstream *file,double val,int precision){
    std::stringstream stream;

    stream << std::fixed << std::setprecision(precision) << val;
    *file  << stream.str() << "\t";
}

void SaveTupleOfDoublesFile(std::ofstream *file, std::valarray<double> tuple, int precision) {
    std::stringstream stream;
    unsigned long i;
    stream << "{";
    for (i = 0; i < tuple.size()-1; ++i) {
        stream << std::fixed << std::setprecision(precision) << tuple[i] << ",";
    }
    stream << std::fixed << std::setprecision(precision) << tuple[i] << "}";
    *file  << stream.str() << "\t";
}

bool is_double(const std::string& s)
{
    std::istringstream iss(s);
    double d;
    char c;
    return iss >> d && !(iss >> c);
}

ParameterFileEntry stringToParameterFileEntry(std::string str_line) {
    std::string name;
    std::vector<std::string> values;
    SplitString(&str_line, &name, &values);
    return ParameterFileEntry(std::move(name), std::move(values));
}

IterableParameterFileEntry stringToIterableParameterFileEntry(std::string str_line) {
    std::string iterate_id;
    std::string name;
    std::vector<std::string> values;
    SplitString(&str_line, &iterate_id, &name, &values);
    return IterableParameterFileEntry(std::move(iterate_id), std::move(name), std::move(values));
}


void checkConsistencyOfIterationParameters(const std::vector<IterableParameterFileEntry>& entries) {
    bool consistent = std::all_of(
        entries.begin(), entries.end(), [entries] (const IterableParameterFileEntry& x) {
            return x.values.size() == entries.front().values.size();
        }
    );
    if (!consistent) {
        throw std::invalid_argument("Iterate 1 entries do not have same length");
    } else if (entries.front().values.size() < 1) {
        throw std::invalid_argument("No values given for iterate entry");
    }
}

void RemoveHashInString(std::vector<std::string> *string){
    std::vector<std::string> new_string{};
    std::string element;
    for (int i{ 0 }; i < string->size(); i++) {
        element = string->at(i);
        if (element.find("#") != std::string::npos) {
            break;
        }
        new_string.push_back(element);
    }
    *string = new_string;
    return;
}
//Struct methods
SynapseExt::SynapseExt(double distToSoma, double lastSpike, double weight): distToSoma{distToSoma}, lastSpike{lastSpike}, weight{weight}
{}

std::valarray<double> SynapseExt::getIndividualSynapticProfile() const {
    std::valarray<double> ret(4);
    ret[0] = this->distToSoma;
    ret[1] = this->theta;
    ret[2] = this->weight;
    ret[3] = this->lastSpike;
    return ret;
};
SynapseExtBranched::SynapseExtBranched(int distanceFromNode, double lastSpike, double weight, int branchId,int branchPositionId) :SynapseExt(static_cast<double>(distanceFromNode), lastSpike, weight), distanceFromNode{ distanceFromNode }, branchId{ branchId }, branchPositionId {branchPositionId}
{}
std::valarray<double> SynapseExtBranched::getIndividualSynapticProfile() const {
    std::valarray<double> ret(3);
    ret[0] = this->distanceFromNode;
    ret[1] = this->weight;
    ret[2] = this->lastSpike;
    return ret;
};

Branch::Branch(int gap, int branchLength, std::vector<int> anteriorBranches, int branchId):spikedSyn(static_cast<size_t>(branchLength/gap), false), synapticGap{gap}, branchLength{branchLength}, anteriorBranches{anteriorBranches}, branchId{branchId}//,branchSynapseIDs(static_cast<size_t>(branchLength/gap), -1)
{
    //std::iota(uniqueSynapsePositionIDs.begin(),uniqueSynapsePositionIDs.end() , branchId*(branchLength/gap));
}

SubRegion::SubRegion(char regionID, std::vector<int> branchesInRegion): regionID{regionID}, branchesInRegion{branchesInRegion}
{

}