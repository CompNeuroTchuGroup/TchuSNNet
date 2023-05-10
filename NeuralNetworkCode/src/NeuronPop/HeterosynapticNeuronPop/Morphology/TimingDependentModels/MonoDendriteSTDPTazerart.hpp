//
// Created by Saif Ahmed on 15.07.21.
//

#ifndef NEURALNETWORK_MONODENDRITESTDPTAZERART_H
#define NEURALNETWORK_MONODENDRITESTDPTAZERART_H

#include "MonoDendriteSTDP.hpp"

class MonoDendriteSTDPTazerart: public MonoDendriteSTDP {

protected:
    void updateLTP(unsigned long synId) override;
    void updateLTD(unsigned long synId) override;

    double gLTP(double deltaT) const override;
    double gLTD(double deltaT) const override;

    double aLTP(double theta) const override;
    double aLTD(double theta) const override;

    double getDistanceEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const override;
    double getTimingEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const override;

    double muLTP{};
    double sigmaLTP{};
    double alpha{};

    double muLTD{};
    double sigmaLTD{};
    double beta{};

public:
    explicit MonoDendriteSTDPTazerart(GlobalSimInfo* info);
    ~MonoDendriteSTDPTazerart() override = default;
    virtual const std::string GetType() override;

    void SaveParameters(std::ofstream * stream, std::string neuronPreId) override;
    void LoadParameters(std::vector<std::string> *input) override;

};


#endif //NEURALNETWORK_MONODENDRITESTDPTAZERART_H
