//
// Created by Saif Ahmed on 02.08.21.
//

#ifndef NEURALNETWORK_MONODENDRITESTDPTAZERARTRELATIVE_H
#define NEURALNETWORK_MONODENDRITESTDPTAZERARTRELATIVE_H

#include "MonoDendriteSTDP.hpp"

class MonoDendriteSTDPTazerartRelative: public MonoDendriteSTDP  {
protected:
    void updateLTP(unsigned long synId) override;
    void updateLTD(unsigned long synId) override;

    double gLTP(double deltaT) const override;
    double gLTD(double deltaT) const override;

    double aLTP(double theta) const override;
    double aLTD(double theta) const override;

    double getDistanceEffects(const SynapseExt* synA, const SynapseExt* synB) const override;
    double getTimingEffects(const SynapseExt* synA, const SynapseExt* synB) const override;

    double muLTP{};
    double sigmaLTP{};
    double alpha{};

    double muLTD{};
    double sigmaLTD{};
    double beta{};

public:
    explicit MonoDendriteSTDPTazerartRelative(GlobalSimInfo* info);
    ~MonoDendriteSTDPTazerartRelative() override = default;
    std::string getType() override;

    void SaveParameters(std::ofstream * stream) override;
    void LoadParameters(std::vector<std::string> *input) override;
};


#endif //NEURALNETWORK_MONODENDRITESTDPTAZERARTRELATIVE_H


