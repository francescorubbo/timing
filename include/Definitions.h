#ifndef DEFINITIONS_H
#define DEFINITIONS_H

enum distribution {gaussian, pseudoRectangular, crabKissingGaussian, crabKissingSquare};
enum smearMode {Off,Z,T,ZT};

const double LIGHTSPEED = 299792458.;
const double PI = 3.141592653589793238463;

typedef CorrInfo pair<double,double>; // (time,eta)

#endif
