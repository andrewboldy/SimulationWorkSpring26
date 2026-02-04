#ifndef MU2E_RANDOM_MOM_ON_SPHERE_HH
#define MU2E_RANDOM_MOM_ON_SPHERE_HH

//----------------------------------------------------------------------------------

// RandomMomOnSphere.hh
//Written by Andrew Boldy
//University of South Carolina
//Spring 2026

//----------------------------------------------------------------------------------
//Helper function that declares the necessary constants, variables, and functions for RandomMomOnSphere.cc
//There are two instances of this:
// a monoenergetic particle being thrown at a given energy, used in ThreeVector RandomMomOnSphere(double magnitude)
// a back to back particle being thrown at a given energy, used in ThreeVector RandomMomOnSphere(double mag1, double mag2, bool B2B, bool DIFFMAG, bool DEBUG)

//----------------------------------------------------------------------------------



//================================================================
// Inclusions
//================================================================

//Standard Inclusions
#include <cmath>
#include <utility>
#include <algorithm>

//CLHEP Inclusions
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/ThreeVector.h"


//Fermilab Inclusions

//CERN ROOT Inclusions

//Personal Inclusions


namespace CLHEP {class HepRandomEngine;} //one line namespace to declare the HepRandomEngine class from CLHEP

namespace mu2e
{
    struct RandomMomOnSphereParams
    {
      double czMin = -1.0;
      double czMax = 1.0;
      double phiMin = 0.0;
      double phiMax = CLHEP::twopi;

      RandomMomOnSphereParams() = default;

      RandomMomOnSphereParams(double cMin, double cMax, double pMin, double pMax)
        : czMin(cMin), czMax(cMax), phiMin(pMin), phiMax(pMax)
        {
          if (czMin > czMax) {std::swap(czMin,czMax);}
          if (phiMin > phiMax) {std::swap(phiMin,phiMax);}
        }
    }; //end parameter struct

    class RandomMomOnSphere
    {
      public :
        //explicit constructor that assumes the full sphere
        explicit RandomMomOnSphere(CLHEP::HepRandomEngine& eng, double czMin = -1., double czMax = 1., double phiMin = 0., double phiMax = CLHEP::twopi);

        //explicit constructor that takes in the parameters in a parameter setup or config file
        explicit RandomMomOnSphere(CLHEP::HepRandomEngine& eng, const RandomMomOnSphereParams& param);

        //~RandomMomOnSphere() noexcept = default; //no special destructor behavior, need to check on this, not comfortable with destructor

        //setters
        void setCZMin(double cMin) {_czMin=cMin;} //setting the czmin value to a given value
        void setCZMax(double cMax) {_czMax=cMax;} //sets the czmax value to a given value
        void setPhiMin(double pMin) {_phiMin=pMin;} //sets the phimin value to a given value
        void setPhiMax(double pMax) {_phiMax=pMax;} //sets the phimax value to a given value


        //getters
        double czMin() const {return _czMin;}
        double czMax() const {return _czMax;}
        double phiMin() const {return _phiMin;}
        double phiMax() const {return _phiMax;}

        //CLHEP::HepRandomEngine& engine() noexcept {return _randFlat.engine();} unsure what noexcept here means, I think it needs to be checked to make sure it works correctly prior
        CLHEP::HepRandomEngine& engine() {return _randFlat.engine();}

        //public member functions for use outside
        CLHEP::Hep3Vector fire(); //creates a random direcitonal unit 3 vector in the unit sphere
        CLHEP::Hep3Vector fire(double mag); //alternative form of fire that takes in one magnitude
        std::pair<CLHEP::Hep3Vector, CLHEP::Hep3Vector> fire_pair(double mag1, double mag2 = 0., bool B2B = true, bool DIFFMAG = false, bool DEBUG = false); //alternative form of fire that takes in two magnitudes

      private :
      //private declaration of the random flat generation
      CLHEP::HepRandomEngine& eng_;
      CLHEP::RandFlat _randFlat;

      //private variable declarations
      double _czMin;
      double _czMax;
      double _phiMin;
      double _phiMax;


    }; //end the RandomMomOnSphere class

} //end namespace mu2e

#endif /* RandomMomOnSphere_hh */
