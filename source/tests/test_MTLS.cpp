#include <UnitTest++/UnitTest++.h>

#include "../src/MTLS.h"
#include "../src/methods.h"

SUITE(MTLS) {
    /* group for tests MTLS */
    class MBTLSFix {
        public:
            MBTLSFix() {
            } // end constructor

            ~MBTLSFix() {
            } // end deconstructor

            template<class DUMMY> double runMinimization(const unsigned int
                    numParameters, double initial, DUMMY* dummy) {
                /* initialize and run minimization, return function value in
                 * minimum */
                Eigen::VectorXd oldParameters =
                    Eigen::VectorXd::Constant(numParameters, initial);
                dummy->calculate(oldParameters);
                Eigen::MatrixXd hessianInverse =
                    Eigen::MatrixXd::Zero(numParameters, numParameters);
                dummy->setHessian(hessianInverse, oldParameters);
                Eigen::VectorXd derivative = dummy->gradient(oldParameters);
                for (unsigned int m = 0; m < 1000; ++m) {
                    /* minimization loop */
                    Eigen::VectorXd searchDirection = - hessianInverse *
                        derivative;
                    Eigen::VectorXd newParameters = oldParameters +
                        searchDirection * Methods::backtrackingLinesearch(1.0,
                                0.5, 0.5, dummy->calculate(oldParameters),
                                oldParameters, searchDirection, derivative,
                                dummy, &DUMMY::calculate);

                    // calculate new values
                    dummy->calculate(newParameters);
                    derivative = dummy->gradient(newParameters);
                   
                    // keep old values
                    oldParameters = newParameters;
                    dummy->setHessian(hessianInverse, newParameters);

                    if (derivative.norm() < 1e-14) {
                        break;
                    } // end if
                } // end form

                return dummy->calculate(oldParameters);
            } // end function runMinimization
    }; // end class MBTLSFix

    TEST_FIXTURE(MBTLSFix, minimizeBacktrackingLinesearch) {
        /* test backtracking linesearch with a simple DFP method for a function
         * f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy() {};
                virtual ~Dummy() {};

                double calculate(const Eigen::VectorXd r) {
                    /* function f(r) = 0.5 * r^2, with r^2 = sum_n x^2_n*/
                    return 0.5*r.squaredNorm();
                } // end function calculate

                Eigen::VectorXd gradient(const Eigen::VectorXd& r) {
                    /* calculate and return derivative */
                    return r;
                } // end function derivative

                void setHessian(Eigen::MatrixXd& H, const Eigen::VectorXd& r) {
                    /* set hessian (dummy which sets to identity for compliance
                     * with general fixture class MBTLSFix) */
                    H = Eigen::MatrixXd::Identity(r.size(), r.size());
                } // end function hessian
        }; // end class Dummy

        // find minimum and make check
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy *dummy = new Dummy();
            double finalValue = runMinimization<Dummy>(n, 5.0, dummy);
            delete dummy;
            CHECK_CLOSE(0, finalValue, 1e-14);
        } // end forn
    } // end TEST_FIXTURE minimizeBacktrackingLinesearch
    
    TEST_FIXTURE(MBTLSFix, minimizeBacktrackingLinesearch2) {
        /* test backtracking linesearch with a simple DFP method for a more
         * complex function f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                double normalFactor, value;

                Dummy(const unsigned int d) {
                    normalFactor = - 1. / pow(sqrt(2. * M_PI), d);
                };

                virtual ~Dummy() {};

                double calculate(const Eigen::VectorXd parameters) {
                    /* function f(r) being the standard Gaussian */
                    value = normalFactor * exp(-0.5 *
                            parameters.squaredNorm());
                    return value;
                } // end function calculate

                Eigen::VectorXd gradient(const Eigen::VectorXd& r) {
                    /* calculate and return derivative */
                    return - r * value;
                } // end function derivative

                void setHessian(Eigen::MatrixXd& H, const Eigen::VectorXd& r) {
                    /* set hessian */
                    for (unsigned int i = 0; i < H.rows(); ++i) {
                        for (unsigned int j = 0; j < H.cols(); ++j) {
                            if (i == j) {
                                H(i,j) = (r(i)*r(i) - 1) * value;
                            } else {
                                H(i,j) = r(i)*r(j) * value;
                            } // end if
                        } // end forj
                    } // end fori
                    H = H.inverse();
                } // end function hessian
        }; // end class Dummy

        // find minimum for n=1,2,3 and run test
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy* dummy = new Dummy(n);
            double finalValue = runMinimization<Dummy>(n, 0.4, dummy);
            delete dummy;
            CHECK_CLOSE(-1. / pow(sqrt(2. * M_PI), n), finalValue, 1e-14);
        } // end fori
    } // end TEST_FIXTURE minimizeBacktrackingLinesearch2

    class MBTLSBFGSFix {
        public:
            MBTLSBFGSFix() {
            } // end constructor

            ~MBTLSBFGSFix() {
            } // end deconstructor

            template<class DUMMY> double runMinimization(const unsigned int
                    numParameters, double initial, DUMMY* dummy) {
                /* initialize and run minimization, return function value in
                 * minimum */
                Eigen::VectorXd oldParameters =
                    Eigen::VectorXd::Constant(numParameters, initial);
                Eigen::MatrixXd hessianInverse =
                    Eigen::MatrixXd::Identity(numParameters, numParameters);
                Eigen::VectorXd derivative = dummy->gradient(oldParameters);
                Eigen::VectorXd oldDerivative = derivative;
                for (unsigned int m = 0; m < 1000; ++m) {
                    /* minimization loop */
                    Eigen::VectorXd searchDirection = - hessianInverse *
                        oldDerivative;
                    Eigen::VectorXd newParameters = oldParameters +
                        searchDirection * Methods::backtrackingLinesearch(0.7,
                                0.5, 0.5, dummy->calculate(oldParameters),
                                oldParameters, searchDirection, derivative,
                                dummy, &DUMMY::calculate);
                    
                    // set derivative and hessian (with BFGS)
                    dummy->calculate(newParameters);
                    derivative = dummy->gradient(newParameters);
                    Methods::BFGSInverse<double>(hessianInverse, newParameters,
                            oldParameters, derivative, oldDerivative);
                   
                    // save old valiues
                    oldParameters = newParameters;
                    oldDerivative = derivative;

                    // simple break condition
                    if (derivative.norm() <= 1e-14) {
                        break;
                    } // end if
                } // end form

                return dummy->calculate(oldParameters);
            } // end function runMinimization
    }; // end class MBTLSBFGSFix

    TEST_FIXTURE(MBTLSBFGSFix, minimizeBacktrackingLinesearchBFGS) {
        /* test backtracking linesearch with a BFGS method for a more complex
         * function f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                double normalFactor, value;

                Dummy(const unsigned int d) {
                    normalFactor = - 1. / pow(sqrt(2. * M_PI), d);
                };

                virtual ~Dummy() {};

                double calculate(const Eigen::VectorXd parameters) {
                    /* function f(r) being the standard Gaussian */
                    value = normalFactor * exp(-0.5 * parameters.squaredNorm());
                    return value;
                } // end function calculate

                Eigen::VectorXd gradient(const Eigen::VectorXd& r) {
                    /* calculate and return derivative */
                    return - r * value;
                } // end function gradient 
        }; // end class Dummy

        // minimize and check
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy* dummy = new Dummy(n);
            double finalValue = runMinimization<Dummy>(n, 0.4, dummy);
            delete dummy;
            CHECK_CLOSE(-1. / pow(sqrt(2. * M_PI), n), finalValue, 1e-14);
        } // end forn
    } // end TEST_FIXTURE minimizeBacktrackingLinesearchBFGS

    class MMTLSFix {
        public:
            MMTLSFix() {
            } // end constructor

            ~MMTLSFix() {
            } // end deconstructor

            template<class DUMMY> double runMinimization(const unsigned int
                    numParameters, double initial, DUMMY* dummy) {
                /* initialize and run minimization, return function value in
                 * minimum */
                Eigen::VectorXd oldParameters =
                    Eigen::VectorXd::Constant(numParameters, initial);
                double oldValue = dummy->calculate(oldParameters);
                Eigen::VectorXd derivative = dummy->gradient();
                Eigen::MatrixXd hessianInverse =
                    Eigen::MatrixXd::Zero(numParameters, numParameters);
                dummy->setHessian(hessianInverse, oldParameters);
                for (unsigned int m = 0; m < 1000; ++m) {
                    /* minimization loop */
                    Eigen::VectorXd searchDirection = - hessianInverse *
                        derivative;
                    Eigen::VectorXd newParameters = oldParameters +
                        searchDirection *
                        MTLS::linesearchMoreThuente<>(searchDirection,
                                oldParameters, oldValue, dummy,
                                &DUMMY::calculate, &DUMMY::gradient);

                    // calculate new values
                    double newValue = dummy->calculate(newParameters);
                    derivative = dummy->gradient();
                   
                    // keep old values
                    oldValue = newValue;
                    oldParameters = newParameters;
                    dummy->setHessian(hessianInverse, newParameters);

                    if (derivative.norm() < 1e-14) {
                        break;
                    } // end if
                } // end form

                return dummy->calculate(oldParameters);
            } // end function runMinimization
    }; // end class MMTLSFix
    
    TEST_FIXTURE(MMTLSFix, minimizeMoreThuenteLinesearch) {
        /* test More Thuente linesearch with a simple DFP method for a function
         * f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy() {};
                virtual ~Dummy() {};

                Eigen::VectorXd derivative;

                double calculate(const Eigen::VectorXd r) {
                    /* function f(r) = 0.5 * r^2, with r^2 = sum_n x^2_n*/
                    derivative = r;
                    return 0.5*r.squaredNorm();
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function derivative

                void setHessian(Eigen::MatrixXd& H, const Eigen::VectorXd& r) {
                    /* set hessian (dummy which sets to identity for compliance
                     * with general fixture class MBTLSFix) */
                    H = Eigen::MatrixXd::Identity(r.rows(), r.rows());
                } // end function hessian
        }; // end class Dummy

        // find minimum and make check
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy *dummy = new Dummy();
            double finalValue = runMinimization<Dummy>(n, 1.0, dummy);
            delete dummy;
            CHECK_CLOSE(0, finalValue, 1e-14);
        } // end forn
    } // end TEST_FIXTURE minimizeMoreThuenteLinesearch
    
    TEST_FIXTURE(MMTLSFix, minimizeMoreThuenteLinesearch2) {
        /* test More-Thuente linesearch with a simple DFP method for a more
         * complex function f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                double normalFactor, value;
                Eigen::VectorXd derivative;

                Dummy(const unsigned int d) {
                    normalFactor = - 1. / pow(sqrt(2. * M_PI), d);
                };

                virtual ~Dummy() {};

                double calculate(const Eigen::VectorXd r) {
                    /* function f(r) being the standard Gaussian */
                    value = normalFactor * exp(-0.5 * r.squaredNorm());
                    derivative = - r * value;
                    return value;
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function derivative

                void setHessian(Eigen::MatrixXd& H, const Eigen::VectorXd& r) {
                    /* set hessian */
                    for (unsigned int i = 0; i < H.rows(); ++i) {
                        for (unsigned int j = 0; j < H.cols(); ++j) {
                            if (i == j) {
                                H(i,j) = (r(i)*r(i) - 1) * value;
                            } else {
                                H(i,j) = r(i)*r(j) * value;
                            } // end if
                        } // end forj
                    } // end fori
                    H = H.inverse();
                } // end function hessian
        }; // end class Dummy

        // find minimum for n=1,2,3 and run test
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy* dummy = new Dummy(n);
            double finalValue = runMinimization<Dummy>(n, 0.5, dummy);
            delete dummy;
            CHECK_CLOSE(-1. / pow(sqrt(2. * M_PI), n), finalValue, 1e-14);
        } // end forn
    } // end TEST_FIXTURE minimizeMoreThuenteLinesearch2

    class MMTLSBFGSFix {
        public:
            MMTLSBFGSFix() {
            } // end constructor

            ~MMTLSBFGSFix() {
            } // end deconstructor

            template<class DUMMY> double runMinimization(const unsigned int
                    numParameters, double initial, DUMMY* dummy) {
                /* initialize and run minimization, return function value in
                 * minimum */
                Eigen::VectorXd oldParameters =
                    Eigen::VectorXd::Constant(numParameters, initial);
                double oldValue = dummy->calculate(oldParameters);
                Eigen::MatrixXd hessianInverse =
                    Eigen::MatrixXd::Identity(numParameters, numParameters);
                Eigen::VectorXd derivative = dummy->gradient();
                Eigen::VectorXd oldDerivative = derivative;
                for (unsigned int m = 0; m < 100; ++m) {
                    /* minimization loop */
                    Eigen::VectorXd searchDirection = - hessianInverse *
                        oldDerivative;
                    Eigen::VectorXd newParameters = oldParameters +
                        searchDirection *
                        MTLS::linesearchMoreThuente<>(searchDirection,
                                oldParameters, oldValue, dummy,
                                &DUMMY::calculate, &DUMMY::gradient);
                    
                    // set derivative and hessian (with BFGS)
                    dummy->calculate(newParameters);
                    derivative = dummy->gradient();
                    Methods::BFGSInverse<double>(hessianInverse, newParameters,
                            oldParameters, derivative, oldDerivative);
                   
                    // save old valiues
                    oldValue = dummy->calculate(oldParameters);
                    oldParameters = newParameters;
                    oldDerivative = derivative;

                    // simple break condition
                    if (derivative.norm() <= 1e-10) {
                        break;
                    } // end if
                } // end form

                return dummy->calculate(oldParameters);
            } // end function runMinimization
    }; // end class MMTLSBFGSFix
    
    TEST_FIXTURE(MMTLSBFGSFix, minimizeMoreThuenteLinesearchBFGS) {
        /* test backtracking linesearch with a simple DFP method for a function
         * f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy() {};
                virtual ~Dummy() {};

                Eigen::VectorXd derivative;

                double calculate(const Eigen::VectorXd r) {
                    /* function f(r) = 0.5 * r^2, with r^2 = sum_n x^2_n*/
                    derivative = r;
                    return 0.5*r.squaredNorm();
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function derivative
        }; // end class Dummy

        // find minimum and make check
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy *dummy = new Dummy();
            double finalValue = runMinimization<Dummy>(n, 1.0, dummy);
            delete dummy;
            CHECK_CLOSE(0, finalValue, 1e-14);
        } // end forn
    } // end TEST_FIXTURE minimizeMoreThuenteLinesearchBFGS

    TEST_FIXTURE(MMTLSBFGSFix, minimizeMoreThuenteLinesearchBFGS1) {
        /* test More-Thuente linesearch with a BFGS method for a more complex
         * function f(r) */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                double normalFactor, value;
                Eigen::VectorXd derivative;

                Dummy(const unsigned int d) {
                    normalFactor = - 1. / pow(sqrt(2. * M_PI), d);
                };

                virtual ~Dummy() {};

                double calculate(const Eigen::VectorXd r) {
                    /* function f(r) being the standard Gaussian */
                    value = normalFactor * exp(-0.5 * r.squaredNorm());
                    derivative = - r * value;
                    return value;
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* return derivative */
                    return derivative;
                } // end function gradient 
        }; // end class Dummy

        // minimize and check
        for (unsigned int n = 1; n <= 3; ++n) {
            Dummy* dummy = new Dummy(n);
            double finalValue = runMinimization<Dummy>(n, 0.4, dummy);
            delete dummy;
            CHECK_CLOSE(-1. / pow(sqrt(2. * M_PI), n), finalValue, 1e-14);
        } // end forn
    } // end TEST_FIXTURE minimizeMoreThuenteLinesearchBFGS

    TEST_FIXTURE(MMTLSBFGSFix, minimizeMoreThuenteLinesearchBFGS2) {
        /* test More-Thuente linesearch with a BFGS method for a more complex
         * function f(r) defined as the Ackleyn function */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Eigen::VectorXd derivative;

                Dummy() {};
                virtual ~Dummy() {};

                double cosSum(const Eigen::VectorXd& r) {
                    double sum = 0;
                    for (unsigned int i = 0; i < r.size(); ++i) {
                        sum += cos(2*M_PI*r(i));
                    } // end fori
                    return sum;
                } // end function cosSum

                double calculate(const Eigen::VectorXd& r) {
                    /* function f(r) being the Ackleyn function */
                    derivative = Eigen::VectorXd::Zero(r.rows());
                    static const double sqrt2 = sqrt(2);
                    static const double eA = exp(1) + 20;
                    double cSum = exp(0.5*cosSum(r));
                    double exprSum = exp(-0.2/sqrt2*r.norm());
                    for (unsigned int i = 0; i < r.size(); ++i) {
                        derivative(i) = 2 * sqrt2 * r(i) * exprSum +
                            M_PI*sin(2*M_PI*r(i)) * cSum;
                    } // end fori
                    return - 20 * exp(-0.2*sqrt(0.5*r.squaredNorm())) - cSum +
                        eA;
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* return derivative */
                    return derivative;
                } // end function gradient 
        }; // end class Dummy

        // calculate final value for test
        Dummy* dummy = new Dummy();
        double finalValue = runMinimization<Dummy>(2, 0.2, dummy);
        delete dummy;

        CHECK_CLOSE(0, finalValue, 1e-14);
    } // end TEST_FIXTURE minimizeMoreThuenteLinesearchBFGS2

    TEST_FIXTURE(MMTLSBFGSFix, minimizeMoreThuenteLinesearchBFGS3) {
        /* test More-Thuente linesearch with a BFGS method for a more complex
         * function f(r) defined as the Rastrigin function */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy() {};
                virtual ~Dummy() {};

                Eigen::VectorXd derivative;
                
                double cosSum(const Eigen::VectorXd& r) {
                    double sum = 0;
                    for (unsigned int i = 0; i < r.size(); ++i) {
                        sum += cos(2*M_PI*r(i));
                    } // end fori
                    return sum;
                } // end function cosSum

                double calculate(const Eigen::VectorXd& r) {
                    /* function f(r) being the Rastrigin function */
                    derivative = Eigen::VectorXd::Zero(r.size());
                    for (unsigned int i = 0; i < r.size(); ++i) {
                        derivative(i) = 2 * (r(i) + 10 * M_PI *
                                sin(2*M_PI*r(i)));
                    } // end fori
                    return 10 * r.size() + r.squaredNorm() - 10 * cosSum(r);
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function gradient 
        }; // end class Dummy

        // create object
        Dummy* dummy = new Dummy();

        // calculate final value for test
        double finalValue = runMinimization<Dummy>(2, 0.2, dummy);

        // free dummy
        delete dummy;

        CHECK_CLOSE(0, finalValue, 1e-14);
    } // end TEST_FIXTURE minimizeMoreThuenteLinesearchBFGS3
    
    class MMTLSonlyFix {
        public:
            MMTLSonlyFix() {
            } // end constructor

            ~MMTLSonlyFix() {
            } // end deconstructor

            template<class DUMMY> double runMinimization(double initial, DUMMY*
                    dummy) {
                /* initialize and run minimization, return function value in
                 * minimum */
                struct Params {
                    /* struct of parameters */
                    double maxIterations = 100;
                    double mu = 0.1;
                    double eta = 0.1;
                    double delta = 4.0;
                    double bisectWidth = 0.66;
                    double bracketTol = 1e-14;
                    double aMin0 = 0.0;
                    double aMax0 = 1000.0;
                } params;

                Eigen::VectorXd oldx = Eigen::VectorXd::Constant(1, initial);
                double oldValue = dummy->calculate(oldx);
                Eigen::VectorXd derivative = dummy->gradient();
                for (unsigned int m = 0; m < 20; ++m) {
                    /* minimization loop */
                    Eigen::VectorXd searchDirection = - derivative;
                    searchDirection /= searchDirection.norm();
                    Eigen::VectorXd newx = oldx + searchDirection *
                        MTLS::linesearchMoreThuente<>(&params, searchDirection,
                                oldx, oldValue, dummy, &DUMMY::calculate,
                                &DUMMY::gradient);

                    // calculate new values
                    double newValue = dummy->calculate(newx);
                    derivative = dummy->gradient();
                   
                    // keep old values
                    oldValue = newValue;
                    oldx = newx;

                    if (fabs(derivative(0)) < 1e-14) {
                        break;
                    } // end if
                } // end form

                return oldValue;
            } // end function runMinimization
    }; // end class MMTLSonlyFix
    
    TEST_FIXTURE(MMTLSonlyFix, minimizeMTLSonly) {
        /* test More-Thuente linesearch on a 1d problem */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy(double mb) {
                    b = mb;
                    derivative = Eigen::VectorXd::Zero(1);
                };
                virtual ~Dummy() {};

                Eigen::VectorXd derivative;
                double b, fmin;

                void setMin() {
                    /* set minimum */
                    fmin = calculate(Eigen::VectorXd::Constant(1,sqrt(b)));
                } // end function setMin

                double calculate(const Eigen::VectorXd a) {
                    /* function f(a)=-a/(a^2 + b) with b a preset constant */
                    derivative(0) = (a(0)*a(0) - b) / pow(a(0)*a(0)+b, 2);
                    return - a(0)/(a(0)*a(0) + b);
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function gradient 
        }; // end class Dummy

        // run tests 
        double b = 2.0;
        Eigen::ArrayXd initials = Eigen::ArrayXd(7);
        initials << 0.0001, 0.001, 0.1, 4.0, 16.0, 100., 1000.;
        for (unsigned int i = 0; i < initials.size(); ++i) {
            Dummy* dummy = new Dummy(b);
            dummy->setMin();
            double exact = dummy->fmin;
            double finalValue = runMinimization<Dummy>(initials(i), dummy);
            delete dummy;
            CHECK_CLOSE(exact, finalValue, 1e-14);
        } // end fori
    } // end TEST_FIXTURE minimizeMTLSonly
    
    TEST_FIXTURE(MMTLSonlyFix, minimizeMTLSonly1) {
        /* test More-Thuente linesearch on a 1d problem */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy(double mb) {
                    b = mb;
                    derivative = Eigen::VectorXd::Zero(1);
                };
                virtual ~Dummy() {};

                Eigen::VectorXd derivative;
                double b, fmin;

                void setMin() {
                    /* set minimum */
                    fmin = calculate(Eigen::VectorXd::Constant(1,8./5-b));
                } // end function setMin

                double calculate(const Eigen::VectorXd a) {
                    /* function f(a)=(a+b)^5 - 2(a+b)^4 with b a preset
                     * constant */
                    double apb = a(0) + b;
                    derivative(0) = apb*apb*apb * (5*apb - 8);
                    return apb*apb*apb*apb*(apb - 2);
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function gradient 
        }; // end class Dummy

        // run tests 
        double b = 0.004;
        Eigen::ArrayXd initials = Eigen::ArrayXd(7);
        initials << 0.001, 0.4, 0.6, 0.8, 1.0, 1.2, 2.0;
        for (unsigned int i = 0; i < initials.size(); ++i) {
            Dummy* dummy = new Dummy(b);
            dummy->setMin();
            double exact = dummy->fmin;
            double finalValue = runMinimization<Dummy>(initials(i), dummy);
            delete dummy;
            CHECK_CLOSE(exact, finalValue, 1e-14);
        } // end fori
    } // end TEST_FIXTURE minimizeMTLSonly1
    
    TEST_FIXTURE(MMTLSonlyFix, minimizeMTLSonly2) {
        /* test More-Thuente linesearch on a 1d problem with a more complex
         * function */

        // define class with a calculation function to send into linesearch
        class Dummy {
            public:
                Dummy(double mb, double ml) {
                    b = mb;
                    l = ml;
                    derivative = Eigen::VectorXd::Zero(1);
                };
                virtual ~Dummy() {};

                Eigen::VectorXd derivative;
                double b, l, fmin;

                void setMin() {
                    /* set minimum */
                    fmin = calculate(Eigen::VectorXd::Constant(1,1));
                } // end function setMin

                double phi0(double a) {
                    if (a <= 1-b) {
                        return 1 - a;
                    } else if (a >= 1+b) {
                        return a - 1;
                    } else {
                        return 0.5 * (1./b * (a-1)*(a-1) + b);
                    } // end ifelse
                } // end function ph0

                double derPhi0(double a) {
                    if (a <= 1-b) {
                        return -1;
                    } else if (a >= 1+b) {
                        return 1;
                    } else {
                        return (a-1)/b;
                    } // end ifelse
                } // end function derPhi

                double calculate(const Eigen::VectorXd a) {
                    /* function f(a)=(a+b)^5 - 2(a+b)^4 with b a preset
                     * constant */
                    double lpi = M_PI * l;
                    double trigFactor = 0.5*lpi*a(0);
                    derivative(0) = derPhi0(a(0)) + (1-b) * cos(trigFactor);
                    return phi0(a(0)) + 2*(1-b)/lpi * sin(trigFactor);
                } // end function calculate

                Eigen::VectorXd gradient() {
                    /* calculate and return derivative */
                    return derivative;
                } // end function gradient 
        }; // end class Dummy

        // run tests 
        double b = 0.01;
        double l = 39;
//         Eigen::ArrayXd initials = Eigen::ArrayXd(8);
//         initials << 0.0, 0.2, 0.6, 1.2, 2.0, 10, 100, 1000;
        Eigen::ArrayXd initials = Eigen::ArrayXd(1);
        initials << 100;
        for (unsigned int i = 0; i < initials.size(); ++i) {
            Dummy* dummy = new Dummy(b, l);
            dummy->setMin();
            double exact = dummy->fmin;
            double finalValue = runMinimization<Dummy>(initials(i), dummy);
            delete dummy;
            CHECK_CLOSE(exact, finalValue, 1e-14);
        } // end fori
    } // end TEST_FIXTURE minimizeMTLSonly2
} // end SUITE MTLS
