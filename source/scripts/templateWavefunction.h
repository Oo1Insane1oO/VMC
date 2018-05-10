#ifndef LCDEF_H
#define LCDEF_H

#include "../basis/clheaderbasis.h"

#include <Eigen/Dense>

class Slater;

class claCL : public clabasis {
    private:
        bool m_interaction;
        
        Eigen::VectorXd m_laplacianSumVec;

        Slater* slater;

    public:
        claCL(Slater*);
        virtual ~claCL();
        
        void setParameters(const Eigen::VectorXd&);
        void initializeParameters();
        void setInteraction(bool);

        double potentialEnergy();
        double kineticEnergy();

        double calculateWavefunction(const unsigned int&, const unsigned int&);
        
        double gradientExpression(const unsigned int&, const int&, const
                unsigned int&);
        const Eigen::VectorXd& laplacianExpression(const unsigned int&, const
                unsigned int&);
        double variationalDerivativeExpression(const unsigned int&, const
                unsigned int&);

        void reSetAll();
        void initializeMatrices();
        void set(const Eigen::MatrixXd& newPositions);
        void update(const Eigen::VectorXd&,  const unsigned int&);
        void reset(const unsigned int&);
        void resetGradient(const unsigned int&);
        void acceptState(const unsigned int&);
        void acceptGradient(const unsigned int&);

    protected:
        std::vector<double(claCL::*)(const unsigned int&, const unsigned int&,
                const double&)> variationalDerivativeFunctionList;
        
};

#endif /* LCDEF_H */
