#include <UnitTest++/UnitTest++.h>

#include "../src/methods.h"

SUITE(Methods) {
    TEST(reMapIndex) {
        /* check that mentioned function gives correct output */
        Eigen::VectorXd values = Eigen::VectorXd::Zero(3);
        values << 0.1, 0.4, 0.8;
        for (unsigned int i = 0; i < values.size(); ++i) {
            CHECK_EQUAL(i, Methods::reMapIndex(values(i), values.size()));
        } // end fori
    } // end TEST probabilityToDimension

    TEST(determinantRatio) {
        /* check that updating determinant ratio works */
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(6,6);
        Eigen::MatrixXd AInv = A.inverse();
        Eigen::MatrixXd B = A;
        B.row(2) = Eigen::VectorXd::Random(6);
        double ratioExact = B.determinant() / A.determinant();
        double ratio = Methods::determinantRatio(B, AInv, 2);
        CHECK_CLOSE(ratioExact, ratio, 1e-14);
    } // end TEST determinantRatio

    TEST(matrixInverse) {
        /* check that updating inverse works */
        Eigen::MatrixXd A = Eigen::MatrixXd::Random(6,6);
        Eigen::MatrixXd AInv = A.inverse();
        Eigen::MatrixXd B = A;
        B.row(2) = Eigen::VectorXd::Random(6);
        Eigen::MatrixXd BInv = B.inverse();
        double ratio = Methods::determinantRatio(B, AInv, 2);
        Eigen::MatrixXd inverse = Eigen::MatrixXd::Zero(6,6);
        Methods::updateMatrixInverse<Eigen::MatrixXd>(A, B, AInv, inverse,
                ratio, 2);
        for (unsigned int i = 0; i < 6; ++i) {
            CHECK_ARRAY_CLOSE(BInv.row(i), inverse.row(i), 6, 1e-14);
        } // end fori
    } // end TEST matrixInverse

    TEST(refSum) {
        /* check that function refSum returns proper sum */
        Eigen::Matrix<int*, Eigen::Dynamic, 1> vec(3);
        Eigen::VectorXi values(3);
        values << 2, 6, 11;
        for (unsigned int i = 0; i < vec.size(); ++i) {
            vec(i) = &(values(i));
        } // end fori
        int sume = values.sum();
        CHECK_EQUAL(sume, Methods::refSum<int>(vec));
    } // end TEST refSum
    
    TEST(refSumMat) {
        /* check that function refSum returns proper sum */
        Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic> mat(3,3);
        Eigen::MatrixXi values(3,3);
        values << 
            2, 6, 11,
            7, 4, 12,
            9, 6, 21;
        for (unsigned int i = 0; i < mat.cols(); ++i) {
            for (unsigned int j = 0; j < mat.cols(); ++j) {
                mat(i,j) = &(values(i,j));
            } // end forj
        } // end forj
        int sume = values.row(2).sum();
        CHECK_EQUAL(sume, Methods::refSum<int>(mat.row(2).segment(0,
                        mat.rows())));
    } // end TEST refSum
} // end SUITE Methods
