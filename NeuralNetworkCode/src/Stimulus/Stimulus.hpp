#ifndef Stimulus_HPP
#define Stimulus_HPP

#include "../GlobalFunctions.hpp"
#include "../NeuronPopSample.hpp"
#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <ranges>

//Struct used in some of the Stimulus classes to coordinate the input current
struct StepStruct {
    TStepInt    startTimeStep{};
    TStepInt    endTimeStep {LONG_MAX};
    std::vector<double> parameterValues;
} ;
/* class Stimulus is a virtual base class for injecting a determined
 * current into each neuron during each time step.
 * - double current(int neuronId, int populationId) returns the current for the
 *   current time step
 * - void timeStepFinished() needs to be called after one time step is finished.
 *   It updates the stimulus object, if neccessary
 * - get_raw_stimulus(int neuronId, int populationId) returns some non-
 *   normalized version of the input current as specified.
 */
class Stimulus {
    
protected:

    GlobalSimInfo*   infoGlobal;
    std::shared_ptr<NeuronPopSample> neurons;
    std::vector<std::vector<double>> signalMatrix;
    int seed;
    std::mt19937 generator;
    bool userSeed {false};

    bool wellDefined{true};
    virtual void SetSignalMatrix() = 0;
    virtual void PostLoadParameters() = 0;

public:

    Stimulus(std::shared_ptr<NeuronPopSample>  neurons,GlobalSimInfo*  infoGlobal);
    virtual ~Stimulus() = default;

    virtual void        Update(std::vector<std::vector<double>>& synaptic_dV) ;
    virtual std::string GetType() const = 0;
   	virtual double GetScaling(PopInt neuronPop) const = 0;

    virtual void SaveParameters(std::ofstream& wParameterStream) const;
    virtual void LoadParameters(const std::vector<FileEntry>& input);

    double GetSignalMatrixPoint(PopInt neuronPop,NeuronInt neuron) const {return signalMatrix.at(neuronPop).at(neuron);}
};
#endif //Stimulus_HPP
