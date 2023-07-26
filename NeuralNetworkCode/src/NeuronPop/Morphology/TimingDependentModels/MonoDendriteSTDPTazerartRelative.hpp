//
// Created by Saif Ahmed on 02.08.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#ifndef NEURALNETWORK_MONODENDRITESTDPTAZERARTRELATIVE_H
#define NEURALNETWORK_MONODENDRITESTDPTAZERARTRELATIVE_H

#include "MonoDendriteSTDP.hpp"

class MonoDendriteSTDPTazerartRelative: public MonoDendriteSTDP  {
    //This class is missing a complete SP!!
protected:
    void UpdateLTP(signed long synId) override;
    void UpdateLTD(signed long synId) override;

    double gLTP(double deltaT) const override;
    double gLTD(double deltaT) const override;

    double aLTP(double theta) const override;
    double aLTD(double theta) const override;


    double muLTP{};
    double sigmaLTP{};
    double alpha{};

    double muLTD{};
    double sigmaLTD{};
    double beta{};

public:
    explicit MonoDendriteSTDPTazerartRelative(GlobalSimInfo* infoGlobal);
    ~MonoDendriteSTDPTazerartRelative() override = default;
    virtual std::string GetType() const override{return IDstringMonoDendriteSTDPTazerartRelative;}

    void SaveParameters(std::ofstream& wParameterStream, std::string neuronIdentificator) const override;
    void LoadParameters(const std::vector<FileEntry>& parameters) override;
    void CheckParameters(const std::vector<FileEntry>& parameters) override;

};


#endif //NEURALNETWORK_MONODENDRITESTDPTAZERARTRELATIVE_H


