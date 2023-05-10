//
// Created by Saif Ahmed on 15.07.21.
//

#ifndef NEURALNETWORK_MONODENDRITESYMMETRICSTDP_H
#define NEURALNETWORK_MONODENDRITESYMMETRICSTDP_H

#include "MonoDendriteSTDP.hpp"

class MonoDendriteSTDPBiWindow: public MonoDendriteSTDP {

protected:
    void updateLTP(unsigned long synId) override;
    void updateLTD(unsigned long synId) override;

    double gLTP(double deltaT) const override;
    double gLTD(double deltaT) const override;

    double aLTP(double theta) const override;
    double aLTD(double theta) const override;

    double getDistanceEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const override;
    double getTimingEffects(const CoopSynapseSpine* synA, const CoopSynapseSpine* synB) const override;

    double tauLTP{};
    double alpha{};

    double tauLTD{};
    double beta{};

public:
    explicit MonoDendriteSTDPBiWindow(GlobalSimInfo* info);
    ~MonoDendriteSTDPBiWindow() override = default;
    virtual const std::string GetType() override;

    void SaveParameters(std::ofstream * stream, std::string neuronPreId) override;
    void LoadParameters(std::vector<std::string> *input) override;
};


#endif //NEURALNETWORK_MONODENDRITESYMMETRICSTDP_H
