#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "vmc.h"
#include "methods.h"
#include "MTLS.h"

template<class T> class VMC;

template<class T>
class Minimizer {
    private:
        VMC<T>* vmc;

        double A, fmax, fmin, a, w, step, tol, maxStepASGD, prevValue, mu, eta,
               eps, exPol, delta, bestEnergy, bestVariance, maxValue,
               oldEnergy, temp, tempStep, beta, astep, annealingFraction;

        struct ParamsMTLS {
            /* struct of default parameters for line search in CG */
            double maxIterations;
            double mu;
            double eta;
            double delta;
            double bisectWidth;
            double bracketTol;
            double aMin0;
            double aMax0;
        } pMTLS; // end struct p1

        unsigned int m_runMax, maxIterationsLinesearch, prevIdx;

        bool hasSetup;

        std::string method, currMethod;

        Eigen::ArrayXd tVec;

        Eigen::VectorXd searchDirection, bestParameters;

        Eigen::MatrixXd hessianInverse, hessian;

        void (Minimizer::*minimizeFunction)();

        void setParamsASGD() {
            /* set specific parameters used in ASGD method */
            A = 60.0;
            fmax = 1.0;
            fmin = -0.5;
            a = 0.9;
            w = 0.5;
            maxStepASGD = 0.1;

            tVec = Eigen::ArrayXd::Constant(vmc->numParameters, A);
        } // end function setParamsASGD

        void setParamsMTLS() {
            /* set parameters used in line search in CG and BFGS method */
            pMTLS.maxIterations = 10;
            pMTLS.mu = 0.001;
            pMTLS.eta = 0.6;
            pMTLS.delta = 4.0;
            pMTLS.bisectWidth = 0.66;
            pMTLS.bracketTol = 1e-14;
            pMTLS.aMin0 = 0.0;
            pMTLS.aMax0 = 1.0;
        } // end function setParamsMTLS

        void setParamsSABFGS() {
            /* set parameters used in SABFGS method */
            beta = 0.6;
        } // end function setParamsSABFGS
        
        void minimizeSD() {
            /* find optimal variational parameters */
            if (!hasSetup) {
                step = 0.001;

                hasSetup = true;
            } // end if

            // find new values
            Methods::steepestDescent<double>(step, vmc->wf->m_parameters,
                    static_cast<Eigen::VectorXd>(vmc->m_newDerivativeParameters
                        / vmc->m_newDerivativeParameters.norm()));
            vmc->sampler();
        } // end function minimizeSD

        void minimizeSABFGS() {
            /* minimize with Stochastic-Adaptive-BFGS method */
            if (!hasSetup) {
                /* allocate for minimization */
                setParamsSABFGS();

                hessianInverse = Eigen::MatrixXd::Identity(vmc->numParameters,
                        vmc->numParameters);
                hessian = Eigen::MatrixXd::Identity(vmc->numParameters,
                        vmc->numParameters);
                searchDirection = Eigen::VectorXd::Zero(vmc->numParameters);

                hasSetup = true;
            } // end if

            static double incrementVal = 1.01;

            auto setSteps = [this]() {
                /* set parameters used in update */
                delta = sqrt(Methods::innerProd(searchDirection, hessian));
                astep = Methods::innerProd(vmc->m_newDerivativeParameters,
                        hessianInverse) / (delta*delta);
                step = 1 / (1./astep + delta);
            }; // end lambda setSteps
            
            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;

            // set search direction (using derivative from previous iteration)
            searchDirection = -hessianInverse * vmc->m_oldDerivativeParameters;

            /* set parameters used in update */
            setSteps();
            
            // sample new values;
            vmc->sampler(vmc->m_oldParameters + step*searchDirection);

            if (Methods::strongWolfeCondition(beta,
                        vmc->m_newDerivativeParameters,
                        vmc->m_oldDerivativeParameters, searchDirection)) {
                /* Make Stochastic-Adaptive step if Strong Wolfe-Conditions are
                 * met */
                searchDirection = -vmc->m_oldDerivativeParameters;
                setSteps();
            } else {
                /* Update Hessian and its inverse with BFGS */
                Methods::BFGSInverse<double>(hessianInverse,
                        vmc->wf->m_parameters, vmc->m_oldParameters,
                        vmc->m_newDerivativeParameters,
                        vmc->m_oldDerivativeParameters);
                Methods::BFGS<double>(hessian, vmc->wf->m_parameters,
                        vmc->m_oldParameters, vmc->m_newDerivativeParameters,
                        vmc->m_oldDerivativeParameters);
            } // end ifelse

            // update parameters
            vmc->setParameters(vmc->m_oldParameters + step*searchDirection);
        } // end function minimizeSABFGS

        void minimizeBFGS() {
            /* find optimal variational parameters */
            if (!hasSetup) {
                /* allocate for minimization */
                setParamsMTLS();

                hessianInverse = Eigen::MatrixXd::Identity(vmc->numParameters,
                        vmc->numParameters);
                searchDirection = Eigen::VectorXd::Zero(vmc->numParameters);

                hasSetup = true;
            } // end if
            
            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;
           
            // set search direction
            searchDirection = - hessianInverse *
                vmc->m_oldDerivativeParameters;
            double sdNorm = searchDirection.norm();
            if (sdNorm >= 1e-10) {
                /* dont normalize if close to minimum for stability */
                searchDirection /= sdNorm;
            } // end if

            // run linesearch
            double s = MTLS::linesearchMoreThuente<>(&pMTLS, searchDirection,
                    vmc->m_oldParameters, vmc->m_accumulativeValues.energy,
                    vmc, static_cast<double(VMC<T>::*)(const Eigen::VectorXd&,
                        const unsigned int)>(&VMC<T>::sampler),
                    &VMC<T>::getNewDerivativeParameters, 0);

            // sample with new parameter and set derivates and update hessian
            // (inverse) with Broyden-Fletcher-Goldfarb-Shanno method
            vmc->setParameters(vmc->m_oldParameters + s*searchDirection);
            vmc->sampler();
            Methods::BFGSInverse<double>(hessianInverse, vmc->wf->m_parameters,
                    vmc->m_oldParameters, vmc->m_newDerivativeParameters,
                    vmc->m_oldDerivativeParameters);
        } // end function minimizeBFGS

        void minimizeCG() {
            /* find optimal variational parameters with the Polak-Ribiere
             * method method */
            if (!hasSetup) {
                /* initialize */
                setParamsMTLS();
                searchDirection = -vmc->m_newDerivativeParameters;

                hasSetup = true;
            } // end if

            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;
           
            // perform linesearch and update parameters
            double s = MTLS::linesearchMoreThuente<>(&pMTLS, searchDirection,
                    vmc->m_oldParameters, vmc->m_accumulativeValues.energy,
                    vmc, static_cast<double(VMC<T>::*)(const Eigen::VectorXd&,
                        const unsigned int)>(&VMC<T>::sampler),
                    &VMC<T>::getNewDerivativeParameters, 0);
            vmc->setParameters(vmc->m_oldParameters + s*searchDirection);

            // evaluate new derivatives and values and update search direction
            // according to Polak-Ribiere method (forcing b to be such that the
            // search direction is always a descent direction)
            vmc->sampler();
            double b =
                Methods::max((vmc->m_newDerivativeParameters.squaredNorm() -
                            vmc->m_newDerivativeParameters.transpose() *
                            vmc->m_oldDerivativeParameters) /
                        vmc->m_oldDerivativeParameters.squaredNorm(), 0.0);
            searchDirection = b * searchDirection -
                vmc->m_newDerivativeParameters;
        } // end function minimizeCG

        void minimizeASGD() {
            /* find optimal variational parameters with Adaptive Stochastic
             * Gradient Descent */
            if (!hasSetup) {
                /* initialize */
                setParamsASGD();

                hasSetup = true;
            } // end if

            // array of values with damping step-function
            Eigen::ArrayXd f = fmin + (fmax - fmin) / (1 - (fmax/fmin) *
                    exp(-(vmc->m_newDerivativeParameters.array() *
                            vmc->m_oldDerivativeParameters.array() / w)));

            // set variables parameters involved in step
            tVec += f;
            for (unsigned int i = 0; i < tVec.size(); ++i) {
                if (tVec(i) < 0) {
                    tVec(i) = 0;
                } // end if
            } // end fori

            // set steps
            Eigen::ArrayXd stepVec = a / (tVec + A);

            // limit to a max step value
            for (unsigned int i = 0; i < tVec.size(); ++i) {
                if (fabs(stepVec(i)) > maxStepASGD) {
                    stepVec(i) *= maxStepASGD / fabs(stepVec(i));
                } // end if
            } // end fori

            // save old values;
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters = vmc->m_newDerivativeParameters;

            // update parameters and derivatives
            vmc->setParameters(vmc->wf->m_parameters - (stepVec *
                        vmc->m_newDerivativeParameters.array() /
                        vmc->m_newDerivativeParameters.norm()).matrix());
            vmc->sampler();
        } // end function minimizeASGD

        void minimizeSIAN() {
            /* minimize with simulated annealing */
            if (!hasSetup) {
                /* initialize */
                vmc->m_oldDerivativeParameters =
                    vmc->m_newDerivativeParameters;
                bestParameters = vmc->wf->m_parameters;
                bestEnergy = vmc->getEnergy();
                bestVariance = (vmc->getEnergySquared() -
                        bestEnergy*bestEnergy) / vmc->m_maxIterations;
                oldEnergy = vmc->getEnergy();
                maxValue = 10.;
                temp = 1.0;
                tempStep = temp/(m_runMax*annealingFraction);

                hasSetup = true;
            } // end if

            // propose new state in a unit gaussian around parameter
            Eigen::ArrayXd newParameters =
                Eigen::ArrayXd::Zero(vmc->wf->m_parameters.size()) .
                unaryExpr([](double val) {
                        static std::normal_distribution<double> nd(val, 0.1);
                        static std::mt19937_64 rng(std::stoi(std::to_string(
                                        std::chrono::
                                        high_resolution_clock::now() .
                                        time_since_epoch() .
                                        count()).substr(10)));
                        return nd(rng);
                    });
            double energy = vmc->sampler(newParameters);

            // use Metropolis-Test to accept/reject state based on the gradient
            // ratio
            double energyDiff = oldEnergy - energy;
            if (vmc->MetropolisTest(exp(energyDiff/temp))) {
                /* update parameters (accept new state) */
                vmc->m_oldParameters = vmc->wf->m_parameters;
                vmc->m_oldDerivativeParameters =
                    vmc->m_newDerivativeParameters;
                double tmpVarience = (vmc->getEnergySquared() - energy*energy)
                    / vmc->m_maxIterations;
                oldEnergy = energy;
                if ((tmpVarience < bestVariance) && (energy < bestEnergy) &&
                        (vmc->getAcceptance() > 0.4)) {
                    /* make choice based on energy variance */
                    bestParameters = vmc->wf->m_parameters;
                    bestEnergy = energy;
                    bestVariance = tmpVarience;
                } // end if
            } else {
                /* reset */
                vmc->setParameters(vmc->m_oldParameters);
                vmc->m_newDerivativeParameters =
                    vmc->m_oldDerivativeParameters;
            } // end if

            temp -= tempStep;
        } // end function minimizeSian

        void setup() {
            vmc->m_oldParameters = vmc->wf->m_parameters;
            vmc->m_oldDerivativeParameters =
                Eigen::VectorXd::Constant(vmc->numParameters, 1000);
            vmc->m_newDerivativeParameters =
                Eigen::VectorXd::Zero(vmc->numParameters);
        } // end function setup

        void showProgress(const unsigned int &index) {
            /* show progress bar for each process */

            // set progress string and buffer
            const std::string progressPosition = Methods::stringPos(vmc->m_rank, 3)
                + "Progress: [";
            std::string progressBuffer;

            // value determining the frequency of progress updates
            static const double divValue = exp(fmod(4,m_runMax));
            if (!static_cast<int>(fmod(index, Methods::divider(index, m_runMax,
                                divValue)))) {
                progressBuffer = progressPosition;
                Methods::printProgressBar(progressBuffer, (float)(((index==m_runMax-1)
                                ? index : index+1)) / m_runMax, 50, currMethod);
            } // end if
        } // end function showProgress

        inline bool updater(const unsigned int& m) {
            /* wrapper for updating minimization methods in function minimize,
             * return false in case break threshold is reached and false if not
             * (return true by default) */

            double testNorm = vmc->m_newDerivativeParameters.norm() /
                Methods::max(1.0, vmc->wf->m_parameters.norm());
//             double testNorm = vmc->m_newDerivativeParameters.norm();
            if (testNorm <= tol) {
                /* break if values are sufficiently convergent */
                return false;
//             } else if ((m>=prevIdx) && ((m%prevIdx)==0) && (((prevValue -
//                                 vmc->m_accumulativeValues.energy) /
//                             vmc->m_accumulativeValues.energy) < tol)) {
//                 /* break if relative improvement is sufficient */
//                 return false;
            } else {
                /* check for switches */
                if (!currMethod.compare("SIAN") && (m > static_cast<unsigned
                            int>(floor(m_runMax*annealingFraction)))) {
                    /* Stop annealing when exhausted and use current best
                     * parameters */
                    vmc->setParameters(bestParameters);
                    if (!method.compare("BFGSADA") || !method.compare("BFGS"))
                    {
                        /* switch to BFGS */
                        currMethod = "BFGS";
                        minimizeFunction = &Minimizer::minimizeBFGS;
                    } else if (!method.compare("SABFGS")) {
                        /* switch to SABFGS */
                        currMethod = "SABFGS";
                        minimizeFunction = &Minimizer::minimizeSABFGS;
                    } else if (!method.compare("CGADA") ||
                            !method.compare("CG")) {
                        /* switch to CG */
                        currMethod = "CG";
                        minimizeFunction = &Minimizer::minimizeCG;
                    } else {
                        /* switch to ASGD */
                        currMethod = "ASGD";
                        minimizeFunction = &Minimizer::minimizeASGD;
                    } // end ifeif
                    hasSetup = false;
                } // end if

                if ((testNorm <= 10*tol) && (!method.compare("BFGSADA")
                            || !method.compare("CGADA"))) {
                    /* change to adaptive stochastic gradient descent method if
                     * sufficiently close to minimum */
                    currMethod = "ASGD";
                    minimizeFunction = &Minimizer::minimizeASGD;
                    hasSetup = false;
                } // end if
            } // end ifelse

            if (!method.compare("SD")) {
                /* reduce step size when close to minimum for steepest descent
                 * */
                if (testNorm <= 100*tol) {
                    step = 0.005;
                } else if (testNorm <= 10*tol) {
                    step = 0.001;
                } // end ifeif
            } // end ifeifelse

            return true;
        } // end function updater
    
    public:
        Minimizer(VMC<T> *vIn, std::string minimizationMethod, unsigned int
                maxMinimizationSamples, double threshold) {
            /* set vmc object */
            vmc = vIn;

            m_runMax = Methods::divider(vmc->m_rank, maxMinimizationSamples,
                    vmc->m_numprocs);

            tol = threshold;

            method = minimizationMethod;
            currMethod = minimizationMethod;
            prevIdx = 10;

            hasSetup = false;

            // always anneal before running BFGS, CG or ASGD
//             currMethod = "SABFGS";
//             minimizeFunction = &Minimizer::minimizeSABFGS;
            currMethod = "SIAN";
            minimizeFunction = &Minimizer::minimizeSIAN;
            annealingFraction = 0.0;

            // set function pointer to specific minimize method
            if (!minimizationMethod.compare("SD")) {
                /* Steepest Descent */
                currMethod = "SD";
                minimizeFunction = &Minimizer::minimizeSD;
            } // end if
            
            if (!(!minimizationMethod.compare("SD") ||
                  !minimizationMethod.compare("CG") ||
                  !minimizationMethod.compare("BFGS") ||
                  !minimizationMethod.compare("CGADA") ||
                  !minimizationMethod.compare("BFGSADA") ||
                  !minimizationMethod.compare("SABFGS") ||
                  !minimizationMethod.compare("ASGD"))) {
                /* throw error if none of the above is given */
                throw std::runtime_error("Please give method, SD, BFGS, ASGD "
                        "BFGSADA, CG or CGADA");
            } // end ifeifeifelse

        } // end constructor
        
        virtual ~Minimizer () {
        } // end deconstructor

        void minimize(bool showBar=true) {
            /* sample and find optimal parameters */
            
            // tell vmc to set derivatives with respect to parameters when
            // sampling
            vmc->m_setVariationalDerivatives = true;
            
            // show progressbar at 0
            if (showBar) {
                showProgress(0);
            } // end if

            // initialize
            setup();
            vmc->sampler();
            prevValue = vmc->m_accumulativeValues.energy;
            std::ofstream outfile;
            outfile.open("test1.txt");
            for (unsigned int m = 0; m < m_runMax; ++m) {
                /* run minimization */

                // show progress
                if (showBar) {
                    showProgress(m);
                } // end if

                if (!currMethod.compare("SIAN")) {
                    outfile << std::setprecision(14) << vmc->getEnergy() << ""
                        " ";
                }

                // check current state of values and update accordingly
                if (!updater(m)) {
                    break;
                } // end if
               
                // update old value
                if ((m>=prevIdx) && ((m%prevIdx)==0)) {
                    prevValue = vmc->m_accumulativeValues.energy;
                } // end if

                // sample and minimize with specified method
                (this->*minimizeFunction)();

//                     vmc->m_newDerivativeParameters.norm() << "      " <<
                std::cout << std::setprecision(14) <<
                    Methods::stringPos(vmc->m_rank+3, vmc->m_numprocs) <<
                    vmc->getAcceptance() << "       " <<
                    vmc->m_newDerivativeParameters.transpose() << "     " <<
                    vmc->wf->m_parameters.transpose() << "   " <<
                    "    " << vmc->m_accumulativeValues.energy << "   " <<
                    (vmc->m_accumulativeValues.energySquared -
                     vmc->m_accumulativeValues.energy *
                     vmc->m_accumulativeValues.energy) / vmc->m_maxIterations
                    << std::endl;
            } // end form

            // show 100% in end
            if (showBar) {
                showProgress(m_runMax-1);
            } // end if
           
            // tell vmc to not set derivatives with respect to parameters when
            // sampling
            vmc->m_setVariationalDerivatives = false;

            if (outfile.is_open()) {
                outfile.close();
            } // end if
        } // end function minimize

        const Eigen::VectorXd& getDerivative() const {
            /* return new derivative vector */
            return vmc->getNewDerivativeParameters();
        } // end function getDerivative

        const Eigen::MatrixXd& getHessianInverse() const {
            /* return new derivative vector */
            return hessianInverse;
        } // end function getDerivative
};

#endif /* MINIMIZER_H */
