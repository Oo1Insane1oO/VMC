#ifdef VARIATIONALHARTREEFOCK

#include "variationalhartreefock.h"
#include "../slater.h"
#include "../hermite/hermite.h"
#include "../methods.h"

#include <boost/math/special_functions/factorials.hpp>

VariationalHartreeFock::VariationalHartreeFock(Slater* sIn) :
    VariationalHartreeFockBasis::VariationalHartreeFockBasis() {
    slater = sIn;
    m_interaction = true;
    
    // add derivatives with respect to variational parameters to list in given
    // order
    variationalDerivativeFunctionList .
        push_back(&VariationalHartreeFock::variationalDerivativeExpression);
} // end constructor

VariationalHartreeFock::~VariationalHartreeFock() {
} // end deconstructor

void VariationalHartreeFock::setHermite3DMatrix(const unsigned int& p) {
    for (unsigned int d = 0; d < slater->m_dim; ++d) {
        for (unsigned int j = 0; j <
                VariationalHartreeFockBasis::Cartesian::getn().size(); ++j) {
            m_hermite3DMatrix(p,d)(j) = H(m_SnewPositions(p,d), j);
        } // end ford
    } // end forj
} // end function setHermite3DMatrix

void VariationalHartreeFock::setParameters(const Eigen::VectorXd& newParameters) {
    /* update parameters */
    /* Set any extra parameters (i.e alpha) here if needed */
    slater->m_parameters = newParameters;

    aw = newParameters(0)*omega;
    sqrtaw = sqrt(aw);

    setHermiteNormalizations();
} // end function setParameters

void VariationalHartreeFock::setInteraction(bool t) {
    /* switch Coulomb interaction on/off */
    m_interaction = t;
} // end function setInteraction

void VariationalHartreeFock::checkIfFullShell() {
    /* make sure shell is filled */
    isFull = false;
    for (unsigned int i = 0; i <
            VariationalHartreeFockBasis::getMagic().size(); ++i) {
        if (slater->m_numParticles==VariationalHartreeFockBasis::getMagic(i)) {
            isFull = true;
            break;
        } // end if
    } // end fori
} // end function checkIfFullShell

void VariationalHartreeFock::initializeParameters(const double& w, const
        unsigned int& L, const Eigen::MatrixXd& coefficientMatrix) {
    /* set number of particles and allocate space for matrices */
    omega = w;
    omegaSq = w * w;
    sqrtOmega = sqrt(w);

    m_numBasis = L;

    // initialize basis (wrapper)
    VariationalHartreeFockBasis::setup(2*L, slater->m_dim);

    m_C = coefficientMatrix.sparseView(1,1e-6);
    m_C.makeCompressed();
   
    // set normalizations
    m_hermiteNormalizations =
        Eigen::ArrayXd::Zero(VariationalHartreeFockBasis::getn().size());
    
    /* fill in any matrices dependant on sizes from basis here */

    // make sure Slater is full-shell (make sure n is a magic number)
    checkIfFullShell();

    // precalculate normalization factors
//     setHermiteNormalizations();
} // end function initializeParameters

double VariationalHartreeFock::variationalDerivativeExpression(const unsigned int& i, const
        unsigned int& j, const double& alphal) {
    /* calculate and return expression in derivative with respect to
     * variational parameter */
    /* fill in */
    double res = 0;
    double lenr = - 0.5 * omega * slater->getNewPosition(i).squaredNorm();
    for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt; ++lIt) {
        double sum = lenr;
        for (unsigned int d = 0; d < slater->m_dim; ++d) {
            const int& n = VariationalHartreeFockBasis::getn(lIt.row(),d);
            if (n==0) {
                continue;
            } // end if
            sum += n/alphal * m_SnewPositions(i,d) *
                m_hermite3DMatrix(i,d)(n-1) / m_hermite3DMatrix(i,d)(n);
        } // end ford
        res += lIt.value() * sum * m_SWavefunctionMatrix(i,lIt.row());
    } // end forlIt
    return res;
} // end function variationalDerivativeExpression

void VariationalHartreeFock::set(const Eigen::MatrixXd& newPositions) {
    /* set function */
    /* set any matrix/vector dependant on the position here. Is will be called
     * in Slater every time corresponding function there is called */
    m_SnewPositions = sqrtaw * newPositions;
    m_SoldPositions = m_SnewPositions;

    for (unsigned int i = 0; i < slater->getNumberOfParticles(); ++i) {
        setHermite3DMatrix(i);
    } // end fori
    m_oldHermite3DMatrix = m_hermite3DMatrix;

    setBasisWavefunction();
    m_SoldWavefunctionMatrix = m_SWavefunctionMatrix;
} // end function set

std::string VariationalHartreeFock::setupDone() {
    /* return message of state of setup */
    std::string possibleN = " ";
    if (isFull) {
        /* return empty message if setup is successfull */
        return "";
    } else {
        /* return message if setup is successfull */
        for (unsigned int i = 0; i <
                VariationalHartreeFockBasis::getMagic().size(); ++i) {
            possibleN +=
                std::to_string(VariationalHartreeFockBasis::getMagic(i)) + " ";
        } // end fori
        return "Slater not full, possible N:" + possibleN;
    } // end ifelse
} // end function setupDone

void VariationalHartreeFock::reSetAll() {
    /* function for reinitializing all matrices except alpha, m_parameters(1),
     * omega and numParticles. */
    /* reinitialize any matrix/vector dependant on the position here. Is will
     * be called in Slater every time corresponding function there is called */
    m_SnewPositions.setZero();
    m_SoldPositions.setZero();
    m_SWavefunctionMatrix.setZero();
    m_SoldWavefunctionMatrix.setZero();
} // end function reSetAll

void VariationalHartreeFock::initializeMatrices() {
    /* initialize matrices with default 0 */
    /* initialize any matrix/vector dependant on the position here. Is will be
     * called in Slater every time corresponding function there is called */
    m_laplacianSumVec = Eigen::VectorXd::Zero(slater->m_numParticles/2);

    m_SnewPositions = Eigen::MatrixXd::Zero(slater->m_numParticles,
            slater->m_dim);
    m_SoldPositions = Eigen::MatrixXd::Zero(slater->m_numParticles,
            slater->m_dim);
    m_SWavefunctionMatrix = Eigen::MatrixXd::Zero(slater->m_numParticles,
            m_numBasis);
    m_SoldWavefunctionMatrix = Eigen::MatrixXd::Zero(slater->m_numParticles,
            m_numBasis);

    m_hermite3DMatrix = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
                      Eigen::Dynamic>::Constant(slater->m_numParticles,
                              slater->m_dim, Eigen::VectorXd::Zero(
                                  VariationalHartreeFockBasis::Cartesian::
                                  getn() . size()));
    m_oldHermite3DMatrix = Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic,
                         Eigen::Dynamic>::Constant(slater->m_numParticles,
                                 slater->m_dim, Eigen::VectorXd::Zero(
                                     VariationalHartreeFockBasis::Cartesian::
                                     getn() . size()));
} // end function setMatricesToZero

void VariationalHartreeFock::setHermiteNormalizations() {
    /* calcualate normalization factors for hermite functions */
    for (unsigned int i = 0; i < m_hermiteNormalizations.size(); ++i) {
        m_hermiteNormalizations(i) = sqrt(sqrtaw / (sqrt(M_PI) * pow(2,i) *
                    boost::math::factorial<double>(i)));
    } // end fori
} // end function setHermiteNormalizations

void VariationalHartreeFock::setBasisWavefunction() {
    /* set m_SWavefunctionMatrix */
    for (unsigned int i = 0; i < slater->m_numParticles; ++i) {
        setBasisWavefunction(i);
    } // end fori
} // end function setBasisWavefunction

void VariationalHartreeFock::setBasisWavefunction(const unsigned int& p) {
    /* set row p in m_SWavefunctionMatrix */
    double expFactor = exp(-0.5*m_SnewPositions.row(p).squaredNorm());
    for (unsigned int l = 0; l < m_numBasis; ++l) {
        m_SWavefunctionMatrix(p,l) = expFactor;
        for (unsigned int d = 0; d < slater->m_dim; ++d) {
            const int& n = VariationalHartreeFockBasis::getn(l,d);
            m_SWavefunctionMatrix(p,l) *= m_hermite3DMatrix(p,d)(n) *
                m_hermiteNormalizations(n);
        } // end ford
    } // end forl
} // end function setBasisWavefunction

void VariationalHartreeFock::update(const Eigen::VectorXd& newPosition, const
        unsigned int& p) {
    /* update position for particle p and calculate new wavefunction, keep old
     * values */
    /* update any matrix/vector dependant on the position or wavefunction here.
     * Is will be called in Slater every time corresponding function there is
     * called */
    m_SoldPositions.row(p) = m_SnewPositions.row(p);
    m_SnewPositions.row(p) = sqrtaw * slater->getNewPosition(p);
    
    m_oldHermite3DMatrix.row(p) = m_hermite3DMatrix.row(p);
    setHermite3DMatrix(p);

    m_SoldWavefunctionMatrix.row(p) = m_SWavefunctionMatrix.row(p);
    setBasisWavefunction(p);
} // end function update

void VariationalHartreeFock::reset(const unsigned int& p) {
    /* set new position and wavefunction to old values for particle p */
    /* reset(revert to old value) any matrix/vector dependant on the position
     * or wavefunction here.  Is will be called in Slater every time
     * corresponding function there is called */
    m_SnewPositions.row(p) = m_SoldPositions.row(p);
    m_SWavefunctionMatrix.row(p) = m_SoldWavefunctionMatrix.row(p);
    
    m_hermite3DMatrix.row(p) = m_oldHermite3DMatrix.row(p);

    m_SWavefunctionMatrix.row(p) = m_SoldWavefunctionMatrix.row(p);
} // end function reset

void VariationalHartreeFock::acceptState(const unsigned int&p) {
    /* accept state and set old positions and matrices accordingly for particle
     * p */
    /* accept(set old to current) any matrix/vector dependant on the the
     * wavefunction. Is will be called in Slater every time corresponding
     * function there is called */
    m_SoldPositions.row(p) =  m_SnewPositions.row(p);
    m_SoldWavefunctionMatrix.row(p) = m_SWavefunctionMatrix.row(p);
    
    m_oldHermite3DMatrix.row(p) = m_hermite3DMatrix.row(p);

    m_SoldWavefunctionMatrix.row(p) = m_SWavefunctionMatrix.row(p);
} // end function acceptState

double VariationalHartreeFock::calculateWavefunction(const unsigned int& p,
        const unsigned int& j) {
    /* calculate and return new wavefunction for particle p in state j */
    /* fill in */
    double res = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt; ++lIt) {
        res += lIt.value() * m_SWavefunctionMatrix(p,lIt.row());
    } // end forl
    return res;
} // end function calculateWavefunction

double VariationalHartreeFock::gradientExpression(const unsigned int& p, const int& j, const
        unsigned int& d) {
    /* calculate gradient expression */
    /* fill in */
    double res = 0.0;
    for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt; ++lIt) {
        const int& nd = VariationalHartreeFockBasis::getn(lIt.row(),d);
        if(nd==0) {
            res -= lIt.value() * m_SnewPositions(p,d) *
                m_SWavefunctionMatrix(p,lIt.row());
        } else {
            res += lIt.value() * (2*nd * m_hermite3DMatrix(p,d)(nd-1) /
                    m_hermite3DMatrix(p,d)(nd) - m_SnewPositions(p,d)) *
                m_SWavefunctionMatrix(p,lIt.row());
        } // end ifelse
    } // end forl
    return res*sqrtaw;
} // end function calculateGradient

const Eigen::VectorXd& VariationalHartreeFock::laplacianExpression(const
        unsigned int& i, const unsigned int& idx) {
    /* calculate and return expression involved in the laplacian */
    /* fill in */
    m_laplacianSumVec.setZero();
    double rlen = m_SnewPositions.row(i).squaredNorm() - slater->m_dim;
    for (unsigned int j = 0; j < idx; ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator lIt(m_C, j); lIt;
                ++lIt) {
            double lsum = rlen;
            for (unsigned int d = 0; d < slater->m_dim; ++d) {
                const int& n = VariationalHartreeFockBasis::getn(lIt.row(),d);
                if (n==0) {
                    continue;
                } else  if (n==1) {
                    lsum -= 4*n/m_hermite3DMatrix(i,d)(n) *
                        m_SnewPositions(i,d)*m_hermite3DMatrix(i,d)(n-1);
                } else {
                    lsum += 4*n/m_hermite3DMatrix(i,d)(n) *
                        ((n-1)*m_hermite3DMatrix(i,d)(n-2) -
                         m_SnewPositions(i,d)*m_hermite3DMatrix(i,d)(n-1));
                } // end if
            } // end ford
            m_laplacianSumVec(j) += lIt.value() * lsum *
                m_SWavefunctionMatrix(i,lIt.row());
        } // end forl
    } // end forj
    m_laplacianSumVec *= aw;
    return m_laplacianSumVec;
} // end function laplacianExpression

double VariationalHartreeFock::potentialEnergy() {
    /* calculate and return potential energy */
    /* fill in */
    double P = 0.5 * omegaSq *
        slater->getNewPositions().rowwise().squaredNorm().sum();
    if (m_interaction) {
        /* run with interaction */
        for (unsigned int i = 0; i < slater->m_numParticles; ++i) {
            for (unsigned int j = i+1; j < slater->m_numParticles; ++j) {
                P += 1. / slater->getNewDistance(i,j);
            } // end forj
        } // end fori
    } // end if
    return P;
} // end function potentialEnergy

#endif
