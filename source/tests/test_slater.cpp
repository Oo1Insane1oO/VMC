#ifdef HARMONICOSCILLATOR

#include <UnitTest++/UnitTest++.h>

#include "../src/slater.h"
#include <Eigen/Dense>
#include <iostream>

SUITE(Quantumdot2D6p) {
    /* group for testing class Quantumdot */
    class QD2D6PFixture {
        public:
            Slater* qd;
            
            unsigned int halfSize;
            Eigen::MatrixXd positions, distances;
            Eigen::VectorXd newPosition, initialParameters;

            QD2D6PFixture() {
                initialParameters = Eigen::VectorXd::Zero(2);
                initialParameters << 0.91, 0.33;
            } // end constructor

            virtual ~QD2D6PFixture() {
            } // end deconstructor

            void setPositions() {
                positions = (-1 +
                        (Eigen::ArrayXXd::Random(qd->getNumberOfParticles(), 2)
                         * 0.5 + 0.5) * 2) * 0.5;
            } // end function setPositions
            
            void setDist() {
                distances = Eigen::MatrixXd::Zero(qd->getNumberOfParticles(),
                        qd->getNumberOfParticles());
                for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
                    for (unsigned int j = i+1; j < qd->getNumberOfParticles();
                            ++j) {
                        distances(i,j) = (positions.row(i) -
                                positions.row(j)).norm();
                        distances(j,i) = distances(i,j);
                    } // end forj
                    for (unsigned int j = 0; j < i;
                            ++j) {
                        distances(i,j) = (positions.row(i) -
                                positions.row(j)).norm();
                        distances(j,i) = distances(i,j);
                    } // end forj
                } // end fori
            } // end function setDistances
    };

    TEST_FIXTURE(QD2D6PFixture, set) {
        /* check that inverse is set properly */
        qd = new Slater(2, 6, initialParameters);
        #ifdef HARMONICOSCILLATOR
            qd->initializeParameters(1.0);
        #endif
        qd->initializeMatrices();
        halfSize = qd->getNumberOfParticles()/2;
        setPositions();
        qd->set(positions);

        for (unsigned int i = 0; i < halfSize; ++i) {
            CHECK_ARRAY_CLOSE(Eigen::MatrixXd::Identity(halfSize,
                        halfSize).row(i),
                    (qd->getWavefunctionMatrix().topRows(halfSize)
                     * qd->getInverseMatrix().topRows(halfSize)).row(i),
                    halfSize, 1e-14);
        } // end fori
        for (unsigned int i = 0; i < halfSize; ++i) {
            CHECK_ARRAY_CLOSE(Eigen::MatrixXd::Identity(halfSize,
                        halfSize).row(i),
                    (qd->getWavefunctionMatrix().bottomRows(halfSize)
                     * qd->getInverseMatrix().bottomRows(halfSize)).row(i),
                    halfSize, 1e-14);
        } // end fori

        delete qd;
    } // end TEST_FIXTURE set

    TEST_FIXTURE(QD2D6PFixture, update2d) {
        /* check that positions are updated properly in 2 dimension */
        qd = new Slater(2, 6, initialParameters);
        #ifdef HARMONICOSCILLATOR
            qd->initializeParameters(1.0);
        #endif
        qd->initializeMatrices();
        halfSize = qd->getNumberOfParticles()/2;
        setPositions();
        qd->set(positions);
        newPosition = (-1 + (Eigen::ArrayXd::Random(2) * 0.5 + 0.5) * 2) * 0.5;
        qd->update(newPosition, 1);
        setDist();
        for (unsigned int j = 2; j < qd->getNumberOfParticles(); ++j) {
            distances(1,j) = (positions.row(1)+newPosition.transpose() -
                    positions.row(j)).norm();
            distances(j,1) = distances(1,j);
        } // end fori
        for (unsigned int j = 0; j < 1; ++j) {
            distances(1,j) = (positions.row(1)+newPosition.transpose() -
                    positions.row(j)).norm();
            distances(j,1) = distances(1,j);
        } // end fori
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            CHECK_ARRAY_CLOSE(distances.row(i), qd->getNewDistance(i),
                    qd->getNumberOfParticles(), 1e-15);
        } // end fori
        CHECK_ARRAY_CLOSE(positions.row(1)+newPosition.transpose(),
                qd->getNewPosition(1), 2, 1e-15);
        CHECK_ARRAY_CLOSE(positions.row(1), qd->getOldPosition(1), 2, 1e-15);
        for (unsigned int i = 0; i < halfSize; ++i) {
            CHECK_ARRAY_CLOSE(Eigen::MatrixXd::Identity(halfSize,
                        halfSize).row(i),
                    (qd->getWavefunctionMatrix().topRows(halfSize)
                     * qd->getInverseMatrix().topRows(halfSize)).row(i),
                    halfSize, 1e-14);
        } // end fori
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            if (i != 1) {
                CHECK_ARRAY_CLOSE(positions.row(i), qd->getNewPosition(i), 2,
                        1e-15);
            } // end if
        } // end fori
        delete qd;
    } // end TEST_FIXTURE update2d

    TEST_FIXTURE(QD2D6PFixture, reset2d) {
        /* check that positions and wavefunction is reset properly */
        qd = new Slater(2, 6, initialParameters);
        #ifdef HARMONICOSCILLATOR
            qd->initializeParameters(1.0);
        #endif
        qd->initializeMatrices();
        halfSize = qd->getNumberOfParticles()/2;
        setPositions();
        qd->set(positions);
        setDist();
        Eigen::MatrixXd inverse = qd->getInverseMatrix();
        double wavefunctionRatio = qd->wavefunctionRatio();
        qd->update(Eigen::VectorXd::Random(2), 0);
        qd->reset(0);
        CHECK_ARRAY_CLOSE(positions.row(0), qd->getNewPosition(0), 2, 1e-15);
        CHECK_ARRAY_CLOSE(positions.row(0), qd->getOldPosition(0), 2, 1e-15);
        CHECK_CLOSE(wavefunctionRatio, qd->wavefunctionRatio(), 1e-15);
        for (unsigned int i = 0; i < halfSize; ++i) {
            CHECK_ARRAY_CLOSE(inverse.topRows(halfSize).row(i),
                    qd->getInverseMatrix().topRows(halfSize).row(i),
                    halfSize, 1e-15);
        } // end fori
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            CHECK_ARRAY_CLOSE(distances.row(i),
                    qd->getNewDistance(i),
                    qd->getNumberOfParticles(), 1e-15);
            CHECK_ARRAY_CLOSE(distances.row(i), qd->getOldDistance(i),
                    qd->getNumberOfParticles(), 1e-15);
        } // end fori
        delete qd;
    } // end TEST_FIXTURE reset2d

    TEST_FIXTURE(QD2D6PFixture, distances) {
        /* check that distance matrix is set correctly */
        qd = new Slater(2, 6, initialParameters);
        #ifdef HARMONICOSCILLATOR
            qd->initializeParameters(1.0);
        #endif
        qd->initializeMatrices();
        halfSize = qd->getNumberOfParticles()/2;
        setPositions();
        qd->set(positions);
        setDist();
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            for (unsigned int j = 0; j < qd->getNumberOfParticles(); ++j) {
                CHECK_CLOSE((qd->getNewPosition(i) -
                            qd->getNewPosition(j)).norm(),
                        qd->getNewDistance(i,j), 1e-15);
            } // end forj
        } // end fori
        const unsigned int p = 4;
        newPosition = (-1 + (Eigen::ArrayXd::Random(2) * 0.5 + 0.5) * 2) * 0.5;
        qd->update(newPosition, p);
        Eigen::MatrixXd locNewPos = positions;
        locNewPos.row(p) += newPosition;
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            for (unsigned int j = 0; j < qd->getNumberOfParticles(); ++j) {
                CHECK_CLOSE((locNewPos.row(i) - locNewPos.row(j)).norm(),
                        qd->getNewDistance(i,j), 1e-15);
            } // end forj
        } // end fori
        qd->reset(p);
        for (unsigned int i = 0; i < qd->getNumberOfParticles(); ++i) {
            for (unsigned int j = 0; j < qd->getNumberOfParticles(); ++j) {
                CHECK_CLOSE((qd->getNewPosition(i) -
                            qd->getNewPosition(j)).norm(),
                        qd->getNewDistance(i,j), 1e-15);
            } // end forj
        } // end fori
        delete qd;
    } // end TEST_FIXTURE distances
} // end SUITE Quantumdot2D6p

#endif
