/*
What is this file for?
*/
#ifndef GUARD_CHI2_FCN_H
#define GUARC_CHI2_FCN_H

#include "Minuit/FCNBase.h"
#include <vector>

class Chi2_Fcn : public FCNBase {
public:
  Chi2_Fcn(const std::vector<double>& meas,
           const std::vector<double>& pos,
           const std::vector<double>& mvar) : theMeasurements(meas),
thePositions(pos),
theMVariances(mvar),
theErrorDef(1.) {}
~GaussFcn() {}
virtual double up() const {return theErrorDef;}
virtual double operator()(const std::vector<double>&) const;
std::vector<double> measurements() const {return theMeasurements;}
std::vector<double> positions() const {return thePositions;}
std::vector<double> variances() const {return theMVariances;}
void setErrorDef(double def) {theErrorDef = def;}
private:
std::vector<double> theMeasurements;
std::vector<double> thePositions;
std::vector<double> theMVariances;
double theErrorDef;
};
#endif //MN_GaussFcn_H_

