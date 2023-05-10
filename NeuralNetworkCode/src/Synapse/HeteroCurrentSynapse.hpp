#ifndef _HETEROCURRENTSYNAPSE_
#define _HETEROCURRENTSYNAPSE_

#include <cassert>

#include "Synapse.hpp"
//#include "../NeuronPop/HeterosynapticNeuronPop/HeteroLIFNeuronPop.hpp"// And removing this include avoids a second include loop
#include "../Connectivity/HeteroDerivedConnectivity/HeteroRandomConnectivity.hpp"
#include "../NeuronPop/HeterosynapticNeuronPop/Morphology/BranchedStructs.hpp"

typedef std::shared_ptr<BaseSynapseSpine> BaseSpinePtr;

class HeteroCurrentSynapse : public Synapse {

protected:

    std::vector<BaseSpinePtr> baseSynapseData{};
    //Branching member variables 
    BranchTargeting synapseTargeting{};

    bool ignoreCouplingStrengthBool{false};
    //void advect_finalize(std::vector<std::vector<double>> * waiting_matrix) override;
    void advectSpikers (std::vector<double>& currents, long spiker) override;
    void advectSpikersSpineCouplingStrength (std::vector<double>& currents, long spiker);
public:

    HeteroCurrentSynapse(NeuronPop* postNeurons, NeuronPop* preNeurons, GlobalSimInfo * info);

    const std::string GetTypeStr() override { return str_heteroSynapse; };
    unsigned long allocateSynapse(unsigned long preId, unsigned long postId) override;

    std::string GetDataHeader(int data_column) override;
	std::string GetUnhashedDataHeader() override;
    std::valarray<double> GetSynapticState(int pre_neuron) override;
    int GetNumberOfDataColumns() override {return 1;};

    void SaveParameters(std::ofstream * stream,std::string id_str) override;
    void LoadParameters(std::vector<std::string> *input) override;

    void advect(std::vector<double> *synaptic_dV) override;
    void advectDefinedCouplingStrength(std::vector<double> *synaptic_dV);
    void advectSpineCouplingStrength(std::vector<double> *synaptic_dV);

    const BranchTargeting& getBranchTarget() const {return synapseTargeting;}// return by const reference for now

    // Testing. Commented out because they are not used in the code
    //friend const std::vector<std::pair<unsigned long, unsigned long>>& getSynapticTargets(HeteroCurrentSynapse&, const unsigned long&);
    //friend std::vector<BaseSynapseSpine> getSynapseData(HeteroCurrentSynapse&);
};

//const std::vector<std::pair<unsigned long, unsigned long>>& getSynapticTargets(HeteroCurrentSynapse&, unsigned long);
//std::vector<BaseSynapseSpine> getSynapseData(HeteroCurrentSynapse&); 

#endif
