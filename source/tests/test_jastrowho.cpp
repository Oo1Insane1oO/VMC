#ifdef HARMONICOSCILLATOR

#include <UnitTest++/UnitTest++.h>
#include <Eigen/Dense>

#include "../src/slaterjastrow.h"

SUITE(Jastrow) {
    /* group for testing class Jastrow */
    class Jastrow2DFix {
        public:
            unsigned int dim;
            unsigned int N;
            Eigen::MatrixXd positions;
            SlaterJastrow* wf;

            Jastrow2DFix() {
                dim = 2;
                N = 20;
                positions = Eigen::MatrixXd::Random(N, dim);
                wf = new SlaterJastrow(dim, N, Eigen::VectorXd::Random(2));
                wf->initializeParameters(1.0);
                wf->initializeMatrices();
            } // end constructor

            ~Jastrow2DFix() {
                delete wf;
            } // end deconstructor
    };
    
    class Jastrow3DFix {
        public:
            unsigned int dim;
            unsigned int N;
            Eigen::MatrixXd positions;
            SlaterJastrow* wf;

            Jastrow3DFix() {
                dim = 3;
                N = 40;
                positions = Eigen::MatrixXd::Random(N, dim);
                wf = new SlaterJastrow(dim, N, Eigen::VectorXd::Random(2));
                wf->initializeParameters(1.0);
                wf->initializeMatrices();
            } // end constructor

            ~Jastrow3DFix() {
                delete wf;
            } // end deconstructor
    };

    TEST_FIXTURE(Jastrow2DFix, checkResetGradient) {
        /* check that gradient is reset properly */
        const unsigned int p = 0;
        wf->set(positions);
        wf->calculateGradient();
        Eigen::MatrixXd originalGradient = wf->getJastrowGradientMatrix();
        wf->update(Eigen::VectorXd::Random(dim), p);
        wf->reset(p);
        for (unsigned int i = 0; i < N; ++i) {
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getJastrowGradient(i), dim, 1e-15);
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getOldJastrowGradient(i), dim, 1e-15);
        } // end fori
    } // end TEST_FIXTURE checkResetGradient
    
    TEST_FIXTURE(Jastrow2DFix, checkAcceptGradient) {
        /* check that gradient is reset properly */
        const unsigned int p = 0;
        wf->set(positions);
        wf->calculateGradient();
        wf->update(Eigen::VectorXd::Random(dim), p);
        wf->calculateGradient(p);
        Eigen::MatrixXd originalGradient = wf->getJastrowGradientMatrix();
        wf->acceptGradient(p);
        for (unsigned int i = 0; i < N; ++i) {
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getJastrowGradient(i), dim, 1e-15);
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getOldJastrowGradient(i), dim, 1e-15);
        } // end fori
    } // end TEST_FIXTURE checkResetGradient

    TEST_FIXTURE(Jastrow3DFix, checkResetGradient) {
        /* check that gradient is reset properly */
        const unsigned int p = 0;
        wf->set(positions);
        wf->calculateGradient();
        Eigen::MatrixXd originalGradient = wf->getJastrowGradientMatrix();
        wf->update(Eigen::VectorXd::Random(dim), p);
        wf->reset(p);
        for (unsigned int i = 0; i < N; ++i) {
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getJastrowGradient(i), dim, 1e-15);
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getOldJastrowGradient(i), dim, 1e-15);
        } // end fori
    } // end TEST_FIXTURE checkResetGradient
    
    TEST_FIXTURE(Jastrow3DFix, checkAcceptGradient) {
        /* check that gradient is reset properly */
        const unsigned int p = 0;
        wf->set(positions);
        wf->calculateGradient();
        wf->update(Eigen::VectorXd::Random(dim), p);
        wf->calculateGradient(p);
        Eigen::MatrixXd originalGradient = wf->getJastrowGradientMatrix();
        wf->acceptGradient(p);
        for (unsigned int i = 0; i < N; ++i) {
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getJastrowGradient(i), dim, 1e-15);
            CHECK_ARRAY_CLOSE(originalGradient.row(i),
                    wf->getOldJastrowGradient(i), dim, 1e-15);
        } // end fori
    } // end TEST_FIXTURE checkResetGradient
} // end SUITE Jastrow

#endif
