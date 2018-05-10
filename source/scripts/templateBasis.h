#ifndef LCDEF_H
#define LCDEF_H

#include <Eigen/Dense>

class claCL {
    private:
        unsigned int m_dim;

    public:
        claCL ();
        virtual ~claCL ();

       void setup();
};

#endif /* LCDEF_H */
