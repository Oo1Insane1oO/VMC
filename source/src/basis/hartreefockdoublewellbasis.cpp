#include "hartreefockdoublewellbasis.h"

HartreeFockDoubleWellBasis::HartreeFockDoubleWellBasis() : DWC(), Cartesian() {
    /* default constructor */
} // end constructor

HartreeFockDoubleWellBasis::HartreeFockDoubleWellBasis(, unsigned int dim) :
    DWC(dim), Cartesian(2*DWC.rows(), dim) {
    /* default constructor */
} // end constructor

HartreeFockDoubleWellBasis::~HartreeFockDoubleWellBasis() {
} // end deconstructor

void HartreeFockDoubleWellBasis::setup(unsigned int dim) {
    /* initiate states */
    DWC::setup(dim);
    Cartesian::setup(2*DWC::rows(), dim);
    Cartesian::restructureStates();
} // end function setup 
