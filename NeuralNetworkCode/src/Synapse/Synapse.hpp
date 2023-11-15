#ifndef SYNAPSE_HPP
#define SYNAPSE_HPP

// class RandomConnectivity;

#include "../GlobalFunctions.hpp"
#include "../Connectivity/Connectivity.hpp"
#include "../Connectivity/RandomConnectivity.hpp"
#include "../Connectivity/PoissonConnectivity.hpp"
#include "../Connectivity/DistanceConnectivity.hpp"
#include "../Connectivity/AdjacencyMatrixConnectivity.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <typeinfo>
#include <cstring>

/*
 * class Synapse is a virtual base class for the population of synapses
 * from one neural population to another
 */
// typedef list
//  typedef void (Synapse::*voidFunctionPtr) (NeuronInt, NeuronInt);
// typedef double (Synapse::*doubleFunctionPtr) (signed long, signed long);

//

class Synapse {
protected:

    std::vector<std::vector<std::pair<NeuronInt, BaseSpinePtr>>> targetSpineList;
    std::vector<std::vector<double>> waitingMatrix; // Refactor
    std::vector<std::vector<int>> DelayDDistribution;       // the list with delays D that are associated with the synapses betweeen the sourceNeuron and its target neurons
    std::vector<std::vector<double>> JCouplingDistribution; // list of synaptic weights J
    // NAMING CONVENTIONS FOR THE CLASS:
    // sourceNeuron == preNeuron, targetNeuron == postNeuron.
    std::mt19937 generator;
    BranchTargeting branchTarget{};
    const PopPtr sourcePop;
    const PopPtr targetPop; // conductance based synapses need access to membrane potential
    const PopInt sourcePopID, targetPopID;

    std::unique_ptr<Connectivity> geometry;
    GlobalSimInfo *infoGlobal;

    int seed{};

    // Connected boolean

    // Connection strength of neurons
    double J{};
    // double      scaledJ{};
    double sigmaJ{};
    double Jpot{};
    double Ppot{0.0};
    // Scaling exponentials
    double scalingFactor{1.0};
    // Synaptic delay
    int Dmin{};
    int Dmax{};
    int delayPerMicrometer{};//Has to be in timesteps per micrometer

    // Cumulated synaptic currents
    double cumulatedDV; // collects all synaptic currents in one time step

    // Connectivity data types:::
 // the list with postsynaptic (or target) neurons and syanpseId (Pair<>) for each neuron of the presynaptic population
    // IMPORTANT: first in pair is target id, second in pair is synapse id, which in ordinary Synapse Objects will be -1.
    double relativeCouplingStrength{1};
                                                            //	bool						HasPot;
    bool HasDdistribution{true};                            // Either (Jpot and Ppot) or Sigma If not HasJDistribution, then all synapses have the same j
    bool HasJDistribution{true};                            // By default, set as true
    bool isRecurrent{false};
    // End of connectivity data
    // HCS data
    bool hasPlasticityModel{false};
    bool ignoreJDParameters{false};
    bool isConnectedBool{false};
    bool userSeed{false}; // by default, seed can be changed

    // End of HCS data
    // Introduce boolean

    void ReadWaitingMatrixEntry(std::vector<double> &synaptic_dV, std::mutex &syndVMutex);

    virtual std::vector<double> AdvectSpikers(NeuronInt spiker) = 0;

    virtual void FillWaitingMatrix(NeuronInt spiker, std::vector<double> &&currents);
    virtual void ResetcumulatedDV() { cumulatedDV = 0.0; };
    // Functions for function pointer GetCurrent (inside Advect spikers)
    double GetCouplingStrength(NeuronInt targetNeuron, NeuronInt sourceNeuron) const;
    double GetDistributionJ(NeuronInt targetNeuronIndex, NeuronInt sourceNeuron) const { return (HasJDistribution) ? JCouplingDistribution.at(sourceNeuron).at(targetNeuronIndex) : J; }
    int GetDistributionD(NeuronInt targetNeuron, NeuronInt sourceNeuron) const { return (HasDdistribution) ? DelayDDistribution.at(sourceNeuron).at(targetNeuron) : Dmax; }
    // double GetCouplingStrengthTrainablePlasticity(signed long sourceNeuron, signed long targetNeuron); //Move to new class, where virtual function goes
    // Necessary mutex locks
    std::mutex _waitingMatrixMutexLock;
    void AllocateSynapseStandard(NeuronInt targetNeuron, NeuronInt sourceNeuron);
    void AllocateSynapseWithPlasticity(NeuronInt targetNeuron, NeuronInt sourceNeuron);

public:
    Synapse(PopPtr targetNeurons, PopPtr sourceNeurons, GlobalSimInfo *infoGlobal);
    virtual ~Synapse() = default;

    // virtual void InitConnect(){}
    //*****************************
    //******* Get Functions *******
    //*****************************
    /**
     * Returns the sum of the coupling strength of an individual neuron
     */
    virtual int GetNoDataColumns() const = 0;
    /**
     * Returns a string that can be used for the data file with the individual
     * parameterName of each column and column number, where data_column is the first
     * column number in the data file.
     */
    virtual std::string GetDataHeader(int dataColumn) = 0;
    virtual std::string GetUnhashedDataHeader() const = 0;
    virtual std::vector<double> GetSynapticState(NeuronInt sourceNeuron) const = 0;

    // std::valarray<double>			GetCurrentcontribution(int pre_neuron);
    double GetRecurrentInput(NeuronInt targetNeuron) const { return waitingMatrix.at((this->infoGlobal->timeStep) % (Dmax + 1)).at(targetNeuron); }

    // double							GetCumulatedDV() const { return cumulatedDV;}
    // int								GetMaxD() const {return D_max;}
    // int								GetMinD() const {return D_min;}
    // double							GetJBase() const {return J;}
    // double							GetJPot() const {return J_pot;}
    // double							GetPPot() const {return P_pot;}

    double GetCumulatedDV() const { return cumulatedDV; }
    int GetMaxD() const { return Dmax; }
    int GetMinD() const { return Dmin; }
    double GetJBase() const { return J; }
    double GetSigmaJ() const { return sigmaJ; }
    double GetJPot() const { return Jpot; }
    double GetPPot() const { return Ppot; }

    virtual std::string GetTypeStr() const = 0;

    std::string GetIdStr() const;
    std::string GetIdStrWithULine() const;

    NeuronInt GetNoSourceNeurons() const { return sourcePop->GetNoNeurons(); }
    NeuronInt GetNoTargetNeurons() const { return targetPop->GetNoNeurons(); }
    NeuronInt GetNoTargetedSynapses(NeuronInt sourceNeuron) const { return static_cast<NeuronInt>(targetSpineList.at(sourceNeuron).size()); }
    PopPtr GetSourceNeuronPop() const { return sourcePop; }
    PopPtr GetTargetNeuronPop() const { return targetPop; }

    // double GetXPositionSourcePop(long neuron) const { return sourcePop->GetXPosition(neuron);}
    // double GetYPositionSourcePop(long neuron) const { return sourcePop->GetYPosition(neuron);}
    // double GetXPositionTargetPop(long neuron) const { return targetPop->GetXPosition(neuron);}
    // double GetYPositionTargetPop(long neuron) const { return targetPop->GetYPosition(neuron);}

    // More extraneous functions
    const BranchTargeting &GetBranchTarget() const { return branchTarget; }

    int GetSeed() const { return seed; }
    PopInt GetSourcePopID() const { return sourcePopID; }
    PopInt GetTargetPopID() const { return targetPopID; }

    bool IsRecurrent() const { return isRecurrent; }
    bool IsConnected() const { return isConnectedBool; }
    //*****************************
    //******* Set Functions *******
    //*****************************
    // virtual void Advect(std::vector<double>& synaptic_dV, std::vector<std::vector<std::vector<double>>> * waiting_matrix);
    // Connectivity functions
    void SetDistributionD();
    void SetDistributionJ(); // For distribution of Js
    virtual void SetSeed(int seed);
    //*****************************
    //******* Misc. Functions *******
    //*****************************

    virtual void ResetWaitingMatrixEntry();
    // virtual void preAdvect();
    void Advect(std::vector<double> &synaptic_dV, std::mutex &_syndVMutex); // NOT SUPPOSED TO BE VIRTUAL
    virtual void LoadParameters(const std::vector<FileEntry> &parameters);
    virtual void SaveParameters(std::ofstream &stream, std::string idString) const;

    virtual void ConnectNeurons();
    // Connectivity functions
    void PostConnectNeurons();

    void WriteConnectivity(const std::string &filename, NeuronInt noNeuronsConnectivity) const;
    void WriteDistributionD(const std::string &filename, NeuronInt noNeuronsWithDelay) const;
    void WriteDistributionJ(const std::string &filename, NeuronInt noNeuronsWithPotentiationJ) const;

    // void Test(){geometry->Test();}
    // Function pointer (look at typedef for type specifications)
    virtual void AllocateSynapse(NeuronInt targetNeuron, NeuronInt sourceNeuron);
    // Functions for if-else tree in allocateSynapse

    bool IsSourceVectorEmpty(NeuronInt sourceNeuron) const { return targetSpineList.at(sourceNeuron).empty(); }
    bool IsTargetLastInVector(NeuronInt targetNeuron, NeuronInt sourceNeuron) const { return targetSpineList.at(sourceNeuron).back().first == targetNeuron; }
};
// std::vector<std::vector<double>> getWaitingMatrix();

#endif // SYNAPSE_HPP
