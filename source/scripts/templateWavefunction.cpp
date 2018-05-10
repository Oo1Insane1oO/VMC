#ifdef LCDEF

#include "clheader.h"
#include "../slater.h"
#include "../methods.h"

claCL::claCL(Slater* sIn) : clabasis::clabasis() {
    slater = sIn;
    m_interaction = true;
    
    // add derivatives with respect to variational parameters to list in given
    // order
    variationalDerivativeFunctionList .
        push_back(&claCL::variationalDerivativeExpression);
} // end constructor

claCL::~claCL() {
} // end deconstructor

void claCL::setParameters(const Eigen::VectorXd& newParameters) {
    /* update parameters */
    /* Set any extra parameters (i.e alpha) here if needed */
    slater->m_parameters = newParameters;
} // end function setParameters

void claCL::setInteraction(bool t) {
    /* switch Coulomb interaction on/off */
    m_interaction = t;
} // end function setInteraction

void claCL::initializeParameters()
    /* set number of particles and allocate space for matrices */

    // initialize basis (wrapper)
    clabasis::setup();
    
    /* fill in any matrices dependant on sizes from basis here */

} // end function initializeParameters

double claCL::variationalDerivativeExpression(const unsigned int& i, const
        unsigned int& j, const double& alphal) {
    /* calculate and return expression in derivative with respect to
     * variational parameter */
    /* fill in */
} // end function variationalDerivativeExpression

void claCL::set(const Eigen::MatrixXd& newPositions) {
    /* set function */
    /* set any matrix/vector dependant on the position here. Is will be called
     * in Slater every time corresponding function there is called */
} // end function set

void claCL::reSetAll() {
    /* function for reinitializing all matrices except alpha, m_parameters(1),
     * omega and numParticles. */
    /* reinitialize any matrix/vector dependant on the position here. Is will
     * be called in Slater every time corresponding function there is called */
} // end function reSetAll

void claCL::initializeMatrices() {
    /* initialize matrices with default 0 */
    /* initialize any matrix/vector dependant on the position here. Is will be
     * called in Slater every time corresponding function there is called */
    m_laplacianSumVec = Eigen::VectorXd::Zero(m_numParticles/2);
} // end function setMatricesToZero

void claCL::update(const Eigen::VectorXd& newPosition, const unsigned int&
        p) {
    /* update position for particle p and calculate new wavefunction, keep old
     * values */
    /* update any matrix/vector dependant on the position or wavefunction here.
     * Is will be called in Slater every time corresponding function there is
     * called */
} // end function update

void claCL::reset(const unsigned int& p) {
    /* set new position and wavefunction to old values for particle p */
    /* reset(revert to old value) any matrix/vector dependant on the position
     * or wavefunction here.  Is will be called in Slater every time
     * corresponding function there is called */
} // end function reset

void claCL::resetGradient(const unsigned int& p) {
    /* set new gradient to old values for particle p */
    /* reset(revert to old value) any matrix/vector dependant on the gradient
     * of the wavefunction. Is will be called in Slater every time
     * corresponding function there is called */
} // end function resetGradient


void claCL::acceptState(const unsigned int&p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    /* accept(set old to current) any matrix/vector dependant on the the
     * wavefunction. Is will be called in Slater every time corresponding
     * function there is called */
} // end function acceptState

void claCL::acceptGradient(const unsigned int&p) {
    /* accept state and set old gradient accordingly for particle
     * p */
    /* accept(set old to current) any matrix/vector dependant on the gradient
     * of the wavefunction. Is will be called in Slater every time
     * corresponding function there is called */
} // end function acceptGradient

double claCL::calculateWavefunction(const unsigned int& p, const unsigned
        int& j) {
    /* calculate and return new wavefunction for particle p in state j */
    /* fill in */
    double res;
    return res;
} // end function calculateWavefunction

double claCL::gradientExpression(const unsigned int& p, const int& j, const
        unsigned int& d) {
    /* calculate gradient expression */
    /* fill in */
} // end function calculateGradient

const Eigen::VectorXd& claCL::laplacianExpression(const unsigned int& i, const
        unsigned int& idx) {
    /* calculate and return expression involved in the laplacian */
    /* fill in */
    m_laplacianSumVec.setZero();
    return m_laplacianSumVec;
} // end function laplacianExpression

double claCL::potentialEnergy() {
    /* calculate and return potential energy */
    /* fill in */
    double P = 0;
    return P;
} // end function potentialEnergy

double claCL::kineticEnergy() {
    /* calculate and return kinetic energy */
    return 0.5 * slater->laplacian();
} // end function kineticEnergy

#endif
