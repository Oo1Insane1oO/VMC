#ifdef HARTREEFOCKDOUBLEWELL

#include "hartreefockdoublewell.h"
#include "../hermite/hermite.h"
#include "../slater.h"
#include "../methods.h"

#include <boost/math/special_functions/factorials.hpp>

HartreeFockDoubleWell::HartreeFockDoubleWell(Slater* sIn) :
    HartreeFockDoubleWellBasis::HartreeFockDoubleWellBasis() {
    slater = sIn;
    m_interaction = true;
} // end constructor

HartreeFockDoubleWell::~HartreeFockDoubleWell() {
} // end deconstructor

void HartreeFock::setHermite3DMatrix(const unsigned int& p) {
    for (unsigned int d = 0; d < slater->m_dim; ++d) {
        for (unsigned int j = 0; j <
                HartreeFockDoubleWellBasis::Cartesian::getn().size(); ++j) {
            m_hermite3DMatrix(p,d)(j) = H(m_SnewPositions(p,d), j);
        } // end ford
    } // end forj
} // end function setHermite3DMatrix

void HartreeFockDoubleWell::setParameters(const Eigen::VectorXd& newParameters) {
    /* update parameters */
    /* Set any extra parameters (i.e alpha) here if needed */
    slater->m_parameters = newParameters;
} // end function setParameters

void HartreeFockDoubleWell::setInteraction(bool t) {
    /* switch Coulomb interaction on/off */
    m_interaction = t;
} // end function setInteraction

void HartreeFockDoubleWell::initializeParameters(const double& w,
        Eigen::MatrixXd coefficientMatrix) {
    /* set number of particles and allocate space for matrices */
    omega = w;
    omegaSq = w * w;
    sqrtOmega = sqrt(w);

    m_numBasis = DWC::C.rows();

    // initialize basis (wrapper)
    HartreeFockDoubleWellBasis::setup(slater->dim);
    
    /* fill in any matrices dependant on sizes from basis here */
    for (unsigned int i = 0; i < coefficientMatrix.rows(); ++i) {
        double value = 0.0;
        for (Eigen::SparseMatrix<double>::InnerIterator itr(DWC::C, i); itr;
                ++itr) {
            itr.valueRef() /= boost::math::factorial<double>(nd) * pos(2,nd) *
                sqrt(M_PI/omega);
            double dval = 1.0;
            for (unsigned int d = 0; d < slater->m_dim; ++d) {
                const int& nd =
                    HartreeFockDoubleWellBasis::Cartesian::getn(lIt.row(), d);
                dval *= 
            } // end ford
            value += itr.value()*itr.value() * dval*dval;
        } // end for itr
        coefficientMatrix.row(i) /= sqrt(value);
    } // end fori
} // end function initializeParameters

double HartreeFockDoubleWell::variationalDerivativeExpression(const unsigned int& i, const
        unsigned int& j, const double& alphal) {
    /* calculate and return expression in derivative with respect to
     * variational parameter */
    /* fill in */
} // end function variationalDerivativeExpression

void HartreeFockDoubleWell::set(const Eigen::MatrixXd& newPositions) {
    /* set function */
    /* set any matrix/vector dependant on the position here. Is will be called
     * in Slater every time corresponding function there is called */
} // end function set

void HartreeFockDoubleWell::reSetAll() {
    /* function for reinitializing all matrices except alpha, m_parameters(1),
     * omega and numParticles. */
    /* reinitialize any matrix/vector dependant on the position here. Is will
     * be called in Slater every time corresponding function there is called */
} // end function reSetAll

void HartreeFockDoubleWell::initializeMatrices() {
    /* initialize matrices with default 0 */
    /* initialize any matrix/vector dependant on the position here. Is will be
     * called in Slater every time corresponding function there is called */
    m_laplacianSumVec = Eigen::VectorXd::Zero(slater->m_numParticles/2);
} // end function setMatricesToZero

void HartreeFockDoubleWell::update(const Eigen::VectorXd& newPosition, const unsigned int&
        p) {
    /* update position for particle p and calculate new wavefunction, keep old
     * values */
    /* update any matrix/vector dependant on the position or wavefunction here.
     * Is will be called in Slater every time corresponding function there is
     * called */
} // end function update

void HartreeFockDoubleWell::reset(const unsigned int& p) {
    /* set new position and wavefunction to old values for particle p */
    /* reset(revert to old value) any matrix/vector dependant on the position
     * or wavefunction here.  Is will be called in Slater every time
     * corresponding function there is called */
} // end function reset

void HartreeFockDoubleWell::resetGradient(const unsigned int& p) {
    /* set new gradient to old values for particle p */
    /* reset(revert to old value) any matrix/vector dependant on the gradient
     * of the wavefunction. Is will be called in Slater every time
     * corresponding function there is called */
} // end function resetGradient


void HartreeFockDoubleWell::acceptState(const unsigned int&p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    /* accept(set old to current) any matrix/vector dependant on the the
     * wavefunction. Is will be called in Slater every time corresponding
     * function there is called */
} // end function acceptState

void HartreeFockDoubleWell::acceptGradient(const unsigned int&p) {
    /* accept state and set old gradient accordingly for particle
     * p */
    /* accept(set old to current) any matrix/vector dependant on the gradient
     * of the wavefunction. Is will be called in Slater every time
     * corresponding function there is called */
} // end function acceptGradient

double HartreeFockDoubleWell::calculateWavefunction(const unsigned int& p, const unsigned
        int& j) {
    /* calculate and return new wavefunction for particle p in state j */
    /* fill in */
    double res;
    return res;
} // end function calculateWavefunction

double HartreeFockDoubleWell::gradientExpression(const unsigned int& p, const int& j, const
        unsigned int& d) {
    /* calculate gradient expression */
    /* fill in */
} // end function calculateGradient

const Eigen::VectorXd& HartreeFockDoubleWell::laplacianExpression(const unsigned int& i, const
        unsigned int& idx) {
    /* calculate and return expression involved in the laplacian */
    /* fill in */
    m_laplacianSumVec.setZero();
    return m_laplacianSumVec;
} // end function laplacianExpression

double HartreeFockDoubleWell::potentialEnergy() {
    /* calculate and return potential energy */
    /* fill in */
    double P = 0;
    return P;
} // end function potentialEnergy

double HartreeFockDoubleWell::kineticEnergy() {
    /* calculate and return kinetic energy */
    return 0.5 * slater->laplacian();
} // end function kineticEnergy

#endif
