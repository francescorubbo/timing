#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <utility>

enum distribution {gaussian, pseudoRectangular, crabKissingGaussian, crabKissingSquare};
enum smearMode {Off,Z,T,ZT};

const double LIGHTSPEED = 299792458.;
const double PI = 3.141592653589793238463;

typedef std::pair<double,double> CorrInfo; // (time,eta)

#endif
