//
//  GlobalFunctions.c
//  NeuralNetworkCode
//
//  Created by Andreas Nold on 13/01/2017.
//  Copyright © 2017 Andreas Nold. All rights reserved.
//
#include "GlobalFunctions.hpp"

DendriticSubRegion::DendriticSubRegion(char regionID, std::vector<int> branchesInRegion): regionID{regionID}, branchesInRegion{branchesInRegion} {
    //This is still under work
}
void FileEntry::RemoveCommentsInValues(char commentCharacter){
    std::vector<std::string> uncommentedValues;
    for (std::string& storedValue: parameterValues){
        if (storedValue.find(commentCharacter) != std::string::npos){
            break;
        } else {
            uncommentedValues.push_back(storedValue);
        }
    }
    this->parameterValues=uncommentedValues;
}
// void MultiplyVector (std::vector<signed long> &vector, signed long value)
// {
//     for(std::vector<signed long>::iterator iterator = vector.begin(); iterator != vector.end(); ++iterator)
//     {
//        iterator =iterator*value;
//     }
// }


// void MultiplyVector (std::vector<double> &vector, double value)
// {
//     for(std::vector<double>::iterator iterator = vector.begin(); iterator != vector.end(); ++iterator)
//     {
//        iterator =iterator*value;
//     }
// }

// void TestWritingFile(std::string filename)
// {
//     // if this test file does not appear in the target directory: stop the
//     // simulation and check the directoryPath.
//     std::ofstream testStream(filename);
//     testStream
//     << "Hello World!\n"
//     << "This is a serious Panda test!\n"
//     << "You passed :-)\n"
//     << "Your directory path was correct.\n"
//     << std::endl;
//     testStream.close();
// }


// void FilterStringVector(std::vector<std::string>& fullString,std::string token,std::vector<std::string>& filteredString)
// {
//     filteredString.clear();
//     for(std::vector<std::string>::iterator iterator = (fullString).begin(); iterator != (fullString).end(); ++iterator) {
//         if(iterator->find(token) != std::string::npos){
//             filteredString.push_back(*iterator);
//         }
//     }
// }

std::vector<FileEntry> FilterStringEntries(const std::vector<FileEntry>& fullEntryVector,std::string filter) {
    std::vector<FileEntry> filteredEntries;
    std::copy_if(fullEntryVector.begin(), fullEntryVector.end(), std::back_inserter(filteredEntries), [filter](const FileEntry& entry){return entry.parameterNameContains(filter);});
    return filteredEntries;
}



std::string getPathToInputFile(std::string& inputFileAddress, bool Windows)  {

    std::istringstream stringStream(inputFileAddress);
    std::string        token;
    std::string        pathTo_inpuFile = "";
    char tokenizerChar = ' ';
    if (Windows)
        tokenizerChar = '\\';
    else
        tokenizerChar = '/';

    while (std::getline(stringStream, token, tokenizerChar)) {
        if (!stringStream.eof()) {
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

// void SplitString(std::string& fullString, std::vector<std::string>& parameterValues) {
// //Converts from string to values (without name involved)
//     std::replace(fullString.begin(), fullString.end(), '\t', ' '); //Replace tabs by spaces, then search for tabs
//     std::replace(fullString.begin(), fullString.end(), '\r', ' ');  //Replace additional end-of-line \r-symbol

//     std::istringstream stringStream(fullString);
//     std::string        token;
//     parameterValues.clear();

//     while (std::getline(stringStream, token, ' ')) {
//         if (token.length() != 0) {
//             parameterValues.push_back(token);
//         }
//     }

// }

std::vector<std::string> SplitStringToValues(std::string fullString) {
//Converts from string to values (without name involved)
    std::vector<std::string> parameterValues;
    std::replace(fullString.begin(), fullString.end(), '\t', ' '); //Replace tabs by spaces, then search for tabs
    std::replace(fullString.begin(), fullString.end(), '\r', ' ');  //Replace additional end-of-line \r-symbol

    std::istringstream stringStream(fullString);
    std::string        token;

    while (std::getline(stringStream, token, ' ')) {
        if (token.length() != 0) {
            parameterValues.push_back(token);
        }
    }
    return parameterValues;
}

FileEntry SplitStringToEntry(std::string fullString)
{
    // Converts from string to values
    std::string parameterName;
    std::vector<std::string> parameterValues;
    std::replace(fullString.begin(),fullString.end(),'\t',' '); //Replace tabs by spaces, then search for tabs
    std::replace(fullString.begin(),fullString.end(),'\r',' ');  //Replace additional end-of-line \r-symbol

    std::istringstream stringStream(fullString);
    std::string        token;

    std::getline(stringStream,parameterName,' '); //put name in name ref
    while(std::getline(stringStream,token,' ')){
        if(token.length() != 0){
            parameterValues.push_back(token);
        }
    }
    return FileEntry{parameterName, parameterValues};
}

IterableFileEntry SplitStringToIterableEntry(std::string fullString) {
    std::string iterateID;
    std::string parameterName;
    std::vector<std::string> parameterValues;
    //Converts from string to values with iterate in the paramline
    std::replace(fullString.begin(),fullString.end(),'\t',' '); //Replace tabs by spaces, then search for tabs
    std::replace(fullString.begin(),fullString.end(),'\r',' ');  //Replace additional end-of-line \r-symbol

    std::istringstream stringStream(fullString);
    std::string        token;

    std::getline(stringStream,iterateID,' '); //put first two entries in their refs
    std::getline(stringStream,parameterName,' ');
    while(std::getline(stringStream,token,' ')){
        if(token.length() != 0){
            parameterValues.push_back(token);
        }
    }
    return IterableFileEntry{iterateID, parameterName, parameterValues};
}

void SaveDoubleFile(std::ofstream& file,double value, int precision) {
    std::stringstream tempStream; //I think the idea of using a stream is to avoid individual writing on the file to slow down the code. As far as I am aware this is already dealt with implicitly in the file stream objects.

    tempStream << std::fixed << std::setprecision(precision) << value;
    file  << tempStream.str() << "\t";
}

// void SaveTupleOfDoublesFile(std::ofstream& file, std::valarray<double> tuple, int precision)  {
//     std::stringstream tempStream; //I think the idea of using a stream is to avoid individual writing on the file to slow down the code. As far as I am aware this is already dealt with implicitly in the file stream objects.
//     size_t dataEntry;
//     tempStream << "{";
//     for (dataEntry = 0; dataEntry < tuple.size()-1; ++dataEntry) {
//         tempStream << std::fixed << std::setprecision(precision) << tuple[dataEntry] << ",";
//     }
//     tempStream << std::fixed << std::setprecision(precision) << tuple[dataEntry] << "}";
//     file  << tempStream.str() << "\t";
// }
void SaveTupleOfDoublesFile(std::ofstream& file, std::vector<double> vector, int precision) {
    std::stringstream tempStream; //I think the idea of using a stream is to avoid individual writing on the file to slow down the code. As far as I am aware this is already dealt with implicitly in the file stream objects.
    size_t dataEntry;
    tempStream << "{"<< std::fixed << std::setprecision(precision);
    for (dataEntry = 0; dataEntry < vector.size()-1; ++dataEntry) {
        tempStream << vector.at(dataEntry) << ",";
    }
    tempStream << vector.at(dataEntry) << "}";
    file  << tempStream.str() << "\t";
}

bool isDouble(const std::string& readString) {
    std::istringstream convertedStringStream(readString);
    double castDouble;
    char castChar;
    return convertedStringStream >> castDouble && !(convertedStringStream >> castChar);
}

size_t IsIterateParamConsistent(FileEntry entry, IterableFileEntry iterateEntry) {
    if (iterateEntry.parameterValues.size()%entry.parameterValues.size()!=0){
        throw "The iterate parameter is not consistent with the introduced parameter.";
    } else {
        return entry.parameterValues.size();
    }
}

signed int MinIterateParameterSize(std::vector<IterableFileEntry> iterateEntries) {
    return static_cast<signed int>(std::min_element(iterateEntries.begin(), iterateEntries.end(), [](const IterableFileEntry& entry1,const IterableFileEntry& entry2){
        return entry1.parameterValues.size()<entry2.parameterValues.size();
    })->parameterValues.size());
}

/*template <typename T>
T stringToFileEntry(std::string readStringLine)
{
    std::string parameterName;
    std::vector<std::string> parameterValues;
    SplitString(&readStringLine, &parameterName, &parameterValues);
    return T(std::move(parameterName), std::move(parameterValues));
}*/


// void CheckConsistencyOfIterationParameters(const std::vector<IterableFileEntry>& entries) {
//     bool consistent = std::all_of( entries.begin(), entries.end(), [entries] (const IterableFileEntry& iterableFileEntry) {
//         return (iterableFileEntry.parameterValues.size() == entries.front().parameterValues.size());
//     });
//     if (entries.front().parameterValues.size() < 1) {
//         throw "No values were given for iterate entry";
//     } else if (!consistent) {
//         throw "Iterate_n entries do not have same length for some n";
//     }
// }

// void RemoveCommentInString(std::vector<std::string>& commentedString, char commentCharacter){
//     std::vector<std::string> uncommentedString{};
//     std::string element;
//     for (int i{ 0 }; i < commentedString.size(); i++) {
//         element = commentedString.at(i);
//         if (element.find(commentCharacter) != std::string::npos) {
//             break;
//         }
//         uncommentedString.push_back(element);
//     }
//     commentedString = uncommentedString;
// }
void threadsafe::put_time(time_t timeObj, const char *formatString, std::stringstream &outputString) {
    std::lock_guard<std::mutex> _lockGuard(_timeMutex);
    outputString<<std::put_time(std::localtime(&timeObj), formatString);
}
