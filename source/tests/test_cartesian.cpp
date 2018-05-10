#include <UnitTest++/UnitTest++.h>

#include "../src/basisfunctions/cartesian.h"
#include <Eigen/Dense>

SUITE(Cartesian3D) {
    /* group for testing class Cartesian in 3 dimensions */
    class QDBasis3DFix {
        public:
            Cartesian* b = new Cartesian();
            QDBasis3DFix() {
                b->setup(40, 3);
            } // end constructor
            ~QDBasis3DFix() {
                delete b;
            } // end deconstructor
    };

    TEST_FIXTURE(QDBasis3DFix, checkStates) {
        /* check that closed-shell is reached for first 4 levels */
        unsigned int numStates = b->getStates().rows();
        CHECK_EQUAL(40, numStates);
    } // end TEST_FIXTURE checkStates

    TEST_FIXTURE(QDBasis3DFix, checkEnergy) {
        /* check that energy is correct */
        Eigen::VectorXi states = b->getE();
        Eigen::VectorXi exact = Eigen::VectorXi(4);
        exact << 1, 2, 3, 4;
        CHECK_ARRAY_EQUAL(exact, states, 4);
    } // end TEST_FIXTURE checkEnergy

    TEST_FIXTURE(QDBasis3DFix, checkLevels) {
        /* check that n-values are correct for 40, 3 particles */
        Eigen::VectorXi n = b->getn();
        Eigen::VectorXi exact = Eigen::VectorXi(4);
        exact << 0, 1, 2, 3;
        CHECK_ARRAY_EQUAL(exact, n, 4);
    } // end TEST_FIXTURE checkLevels

    TEST_FIXTURE(QDBasis3DFix, checkMagic) {
        /* check that magic numbers are correct */
        Eigen::VectorXi M = b->getMagic();
        Eigen::VectorXi exact = Eigen::VectorXi(4);
        exact << 2, 8, 20, 40;
        CHECK_ARRAY_EQUAL(exact, M, 4);
    } // end TEST_FIXTURE checkMagic
    
    TEST_FIXTURE(QDBasis3DFix, checkRestructure) {
        /* check that spin states are in increasing order after restructuring
         * */
        b->restructureStates();
        const unsigned int half = b->getStates().rows()/2;
        for (unsigned int i = 0; i < half; ++i) {
            CHECK_EQUAL(-1, *(b->getStates(i)(4)));
        } // end fori
        for (unsigned int i = half; i < b->getStates().rows(); ++i) {
            CHECK_EQUAL(1, *(b->getStates(i)(4)));
        } // end fori
    } // end TEST_FIXTURE checkRestructure
} // end SUITE quantumdotBasis3D

SUITE(Cartesian2D) {
    /* group for testing class quantumdotBasis in 2 dimensions */
    class QDBasis2DFix {
        public:
            Cartesian* b = new Cartesian();
            QDBasis2DFix() {
                b->setup(42, 2);
            } // end constructor
            ~QDBasis2DFix() {
                delete b;
            } // end deconstructor
    };

    TEST_FIXTURE(QDBasis2DFix, checkStates) {
        /* check that closed-shell is reached for first 4 levels */
        unsigned int numStates = b->getStates().rows();
        CHECK_EQUAL(42, numStates);
    } // end TEST_FIXTURE checkStates

    TEST_FIXTURE(QDBasis2DFix, checkEnergy) {
        /* check that energy is correct */
        Eigen::VectorXi states = b->getE();
        Eigen::VectorXi exact = Eigen::VectorXi(6);
        exact << 1, 2, 3, 4, 5, 6;
        CHECK_ARRAY_EQUAL(exact, states, 3);
    } // end TEST_FIXTURE checkEnergy

    TEST_FIXTURE(QDBasis2DFix, checkLevels) {
        /* check that n-values are correct for 20 particles */
        Eigen::VectorXi n = b->getn();
        Eigen::VectorXi exact = Eigen::VectorXi(6);
        exact << 0, 1, 2, 3, 4, 5;
        CHECK_ARRAY_EQUAL(exact, n, 4);
    } // end TEST_FIXTURE checkLevels

    TEST_FIXTURE(QDBasis2DFix, checkMagic) {
        /* check that magic numbers are correct */
        Eigen::VectorXi M = b->getMagic();
        Eigen::VectorXi exact = Eigen::VectorXi(5);
        exact << 2, 6, 12, 20, 42;
        CHECK_ARRAY_EQUAL(exact, M, 4);
    } // end TEST_FIXTURE checkMagic

    TEST_FIXTURE(QDBasis2DFix, checkSumn) {
        /* check that summing n-values is working */
        Eigen::Matrix<int*, Eigen::Dynamic, 1> nvals =
            b->getStates().row(6).segment(0,2);
        int sumn = Methods::refSum(nvals);
        CHECK_EQUAL(sumn, b->getSumn(6));
    } // endf TEST_FIXTURE checkSumn

    TEST_FIXTURE(QDBasis2DFix, checkRestructure) {
        /* check that spin states are in increasing order after restructuring
         * */
        b->restructureStates();
        const unsigned int half = b->getStates().rows()/2;
        for (unsigned int i = 0; i < half; ++i) {
            CHECK_EQUAL(-1, *(b->getStates(i)(3)));
        } // end fori
        for (unsigned int i = half; i < b->getStates().rows(); ++i) {
            CHECK_EQUAL(1, *(b->getStates(i)(3)));
        } // end fori
    } // end TEST_FIXTURE checkRestructure
} // end SUITE quantumdotBasis2D
