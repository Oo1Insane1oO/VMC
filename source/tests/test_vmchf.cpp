#ifdef HARTREEFOCK

#include <UnitTest++/UnitTest++.h>
#include <yaml-cpp/yaml.h>
#include <string>
#include <unistd.h>

#include "../src/examples/run.h"

SUITE(VMCHF) {
    /* group for testing VMC runs with Hartree-Fock basis */
    class fix1 {
        public:
            YAML::Node node;

            fix1() {
            }
            virtual ~fix1() {}

            void setNodes(bool i, double s) {
                node["test"] = true;
                node["jastrow"] = false;
                node["progress"] = false;
                node["importance"] = i;
                node["stepmc"] = s;
            }
    };

    TEST_FIXTURE(fix1, energiesw1N2D2L30) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        node = YAML::LoadFile("tests/inputs/w_1.0_D_2_N_2_L_30.yaml");
        setNodes(false, 1.3);
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30
    
    TEST_FIXTURE(fix1, energiesw1N2D2L30Imp) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        node = YAML::LoadFile("tests/inputs/w_1.0_D_2_N_2_L_30.yaml");
        setNodes(true, 0.01);
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30Imp
    
    TEST_FIXTURE(fix1, energiesw028N2D2L30) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        node = YAML::LoadFile("tests/inputs/w_0.28_D_2_N_2_L_20.yaml");
        setNodes(false, 1.3);
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30
       
    TEST_FIXTURE(fix1, energiesw028N2D2L30Imp) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        node = YAML::LoadFile("tests/inputs/w_0.28_D_2_N_2_L_20.yaml");
        setNodes(true, 2.5);
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30Imp
    
    TEST_FIXTURE(fix1, energiesw1N6D2L42) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        node = YAML::LoadFile("tests/inputs/w1.0_N6_D2_L42.yaml");
        setNodes(false, 1.2);
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30
       
    TEST_FIXTURE(fix1, energiesw1N6D2L42Imp) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        node = YAML::LoadFile("tests/inputs/w1.0_N6_D2_L42.yaml");
        setNodes(true, 0.02);
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30Imp
} // end SUITE VMCHF

#endif
