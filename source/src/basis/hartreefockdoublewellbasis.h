#ifndef HARTREEFOCKDOUBLEWELLBASIS_H
#define HARTREEFOCKDOUBLEWELLBASIS_H

#include "../basisfunctions/cartesian.h"
#include "../basisfunctions/dwc.h"

#include <Eigen/Dense>

class HartreeFockDoubleWellBasis : public DWC, Cartesian {
    private:
        unsigned int m_dim;

    public:
        HartreeFockDoubleWellBasis ();
        HartreeFockDoubleWellBasis (unsigned int, unsigned int);
        virtual ~HartreeFockDoubleWellBasis ();

        void setup(unsigned int, unsigned int);
};

#endif /* HARTREEFOCKDOUBLEWELLBASIS_H */
