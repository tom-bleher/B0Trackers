#ifndef VARIABLES_H
#define VARIABLES_H

#include <array>
#include <vector>
#include <math.h>
#include <cstdint>

struct HitClass{
    double xR;
    double yR;
    double zR;
    
    double xT;
    double yT;
    double zT;
    
    double px;
    double py;
    double pz;
    double p;
    double theta;
    double phi;
 
    int prime;   
    int plane;
    int mod;
    int side;
    uint64_t cellID;
    double eDep;
    double time;
    double pathL;

    double xDet;
    double yDet;
    double zDet;
};

struct TrackClass{
    
    uint64_t cellID;
    double e;
    double charge;
    double x0;
    double y0;
    double slopeX;
    double slopeY;
    double theta;
    double phi;
    double chi2;
    double prime;
};

struct Track {
    std::vector<const HitClass*> hits;
    double chi2 = 0.0;  // can be filled later by a fit function
    // Fit result quantities
    double x0 = 0.0;      // x at z mean
    double y0 = 0.0;      // y at z mean
    double slopeX = 0.0;  // dx/dz
    double slopeY = 0.0;  // dy/dz
    double theta = 0.0;   // polar angle
    double phi   = 0.0;   // azimuthal angle
    int isPrimary = 0;
    int side = -1; // top or bottom Top - 1, bottom - 0
    std::vector<std::array<double, 3>> recoPositions;  // (xR, yR, zR)
    std::vector<std::array<double, 3>> truePositions;  // (xT, yT, zT)
    int charge = 0;
};


struct VertexDataTrue{

	std::vector<int> status;
	std::vector<double> E;
	std::vector<double> vx;
	std::vector<double> vy;
	std::vector<double> vz;
	std::vector<double> px;
	std::vector<double> py;
	std::vector<double> pz;
	std::vector<double> p;
	std::vector<double> theta;
	std::vector<double> phi;
	std::vector<int> pdg;

};

struct Node {
    int planeIndex;
    std::vector<const HitClass*> seed;
};

#endif

