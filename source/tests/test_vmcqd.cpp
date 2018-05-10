#ifdef HARMONICOSCILLATOR

#include <UnitTest++/UnitTest++.h>

#include "../src/slaterjastrow.h"
#include "../src/slater.h"
#include "../src/bruteforce.h"
#include "../src/importanceSampling.h"

SUITE(VMCQD) {
    /* group for testing VMC runs with quantumdot */
    TEST(distances) {
        /* check that distance matrices are updated correctly */
        Eigen::VectorXd initialParameters(2);
        initialParameters << 1, 1e16;
        Slater* qd = new Slater(2,6,initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        Bruteforce<Slater>* vmc = new Bruteforce<Slater>(qd, 1.2,
                initialParameters, 2, 0, 1);

        vmc->sampler();
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            for (unsigned int j = 0; j < qd->getNumberOfParticles(); ++j) {
                CHECK_CLOSE((qd->getNewPosition(i) -
                            qd->getNewPosition(j)).norm(),
                        qd->getNewDistance(i,j), 1e-15);
            } // end forj
        } // end fori
        delete vmc;
        delete qd;
    } // end TEST distances
    
    TEST(distances2) {
        /* check that distance matrices are updated correctly */
        Eigen::VectorXd initialParameters(2);
        initialParameters << 1, 1e16;
        Slater* qd = new Slater(2,6,initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<Slater>* vmc = new ImportanceSampling<Slater>(qd,
                0.001, initialParameters, 2, 0, 1);

        vmc->sampler();
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            for (unsigned int j = 0; j < qd->getNumberOfParticles(); ++j) {
                CHECK_CLOSE((qd->getNewPosition(i) -
                            qd->getNewPosition(j)).norm(),
                        qd->getNewDistance(i,j), 1e-15);
            } // end forj
        } // end fori
        delete vmc;
        delete qd;
    } // end TEST distances

    TEST(unperturbedEnergies2) {
        /* check that energy in case without interaction is correct for 2
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,2,initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        Bruteforce<Slater>* vmc = new Bruteforce<Slater>(qd, 1.2,
                initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(2, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies2

    TEST(unperturbedEnergies6) {
        /* check that energies in case without interaction is correct for 6
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,6, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        Bruteforce<Slater>* vmc = new Bruteforce<Slater>(qd, 1.2,
                initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(10, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies6
    
    TEST(unperturbedEnergies12) {
        /* check that energies in case without interaction is correct for 12
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,12, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        Bruteforce<Slater>* vmc = new Bruteforce<Slater>(qd, 1.2,
                initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(28, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies12
    
    TEST(unperturbedEnergies20) {
        /* check that energies in case without interaction is correct for 20
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,20, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        Bruteforce<Slater>* vmc = new Bruteforce<Slater>(qd, 1.2,
                initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(60, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies20
    
    TEST(unperturbedEnergies2Imp) {
        /* check that energies in case without interaction is correct for 2
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,2, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<Slater>* vmc = new ImportanceSampling<Slater>(qd,
                0.001, initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(2, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies2
    
    TEST(unperturbedEnergies6Imp) {
        /* check that energies in case without interaction is correct for 6
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,6, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<Slater>* vmc = new ImportanceSampling<Slater>(qd,
                0.001, initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(10, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies6
    
    TEST(unperturbedEnergies12Imp) {
        /* check that energies in case without interaction is correct for 12
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,12, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<Slater>* vmc = new ImportanceSampling<Slater>(qd,
                0.01, initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(28, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies12
    
    TEST(unperturbedEnergies20Imp) {
        /* check that energies in case without interaction is correct for 20
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        Slater* qd = new Slater(2,20, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<Slater>* vmc = new ImportanceSampling<Slater>(qd,
                0.01, initialParameters, 1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(60, vmc->getEnergy(), 1e-10);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies20

    TEST(unperturbedEnergies2ImpJast) {
        /* check that energies in case without interaction is correct for 2
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,2, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.001, initialParameters,
                    1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(2, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies2
    
    TEST(unperturbedEnergies6ImpJast) {
        /* check that energies in case without interaction is correct for 6
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,6, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.001, initialParameters,
                    1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(10, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies6
    
    TEST(unperturbedEnergies12ImpJast) {
        /* check that energies in case without interaction is correct for 12
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,12, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.01, initialParameters,
                    1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(28, vmc->getEnergy(), 1e-12);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies12
    
    TEST(unperturbedEnergies20ImpJast) {
        /* check that energies in case without interaction is correct for 20
         * particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 1,
                1e16).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,20, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        qd->setInteraction(false);
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.01, initialParameters,
                    1e4, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(60, vmc->getEnergy(), 1e-10);
        delete vmc;
        delete qd;
    } // end TEST unperturbedEnergies20

    TEST(interactionEnergies2) {
        /* check that energies with interaction is correct for 2 particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 0.99,
                0.40).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,2, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        Bruteforce<SlaterJastrow>* vmc = new Bruteforce<SlaterJastrow>(qd, 1.3,
                initialParameters, 1e5, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(3.0, vmc->getEnergy(), 1e-2);
        delete vmc;
        delete qd;
    } // end TEST interactionEnergies2
    
    TEST(interactionEnergies2Imp) {
        /* check that energies with interaction is correct for 2 particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 0.99,
                0.40).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,2, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.001, initialParameters,
                    1e5, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(3.0, vmc->getEnergy(), 1e-2);
        delete vmc;
        delete qd;
    } // end TEST interactionEnergies2
    
    TEST(interactionEnergies6) {
        /* check that energies with interaction is correct for 2 particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 0.99,
                0.47).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,6, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        Bruteforce<SlaterJastrow>* vmc = new Bruteforce<SlaterJastrow>(qd, 1.3,
                initialParameters, 1e5, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(20.5, vmc->getEnergy(), 0.5);
        delete vmc;
        delete qd;
    } // end TEST interactionEnergies6
    
    TEST(interactionEnergies6Imp) {
        /* check that energies with interaction is correct for 2 particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 0.99,
                0.47).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,6, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.001, initialParameters,
                    1e5, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(20.5, vmc->getEnergy(), 0.5);
        delete vmc;
        delete qd;
    } // end TEST interactionEnergies6Imp
    
    TEST(interactionEnergies12) {
        /* check that energies with interaction is correct for 2 particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 0.68,
                1.06).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,12, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        Bruteforce<SlaterJastrow>* vmc = new Bruteforce<SlaterJastrow>(qd, 1.5,
                initialParameters, 1e5, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(67.1, vmc->getEnergy(), 0.5);
        delete vmc;
        delete qd;
    } // end TEST interactionEnergies6
    
    TEST(interactionEnergies12Imp) {
        /* check that energies with interaction is correct for 2 particles */
        Eigen::VectorXd initialParameters = (Eigen::VectorXd(2) << 0.68,
                1.06).finished();
        SlaterJastrow* qd = new SlaterJastrow(2,12, initialParameters);
        qd->initializeParameters(1.0);
        qd->initializeMatrices();
        ImportanceSampling<SlaterJastrow>* vmc = new
            ImportanceSampling<SlaterJastrow>(qd, 0.01, initialParameters,
                    1e5, 0, 1);

        vmc->sampler();
        CHECK_CLOSE(67.1, vmc->getEnergy(), 0.5);
        delete vmc;
        delete qd;
    } // end TEST interactionEnergies6Imp
} // end SUITE VMCQD

#endif
