#ifdef HARTREEFOCK

#include <UnitTest++/UnitTest++.h>
#include <yaml-cpp/yaml.h>
#include <string>
#include <unistd.h>

#include "../src/examples/run.h"

SUITE(VMCHF) {
    /* group for testing VMC runs with Hartree-Fock basis */
    TEST(energiesw1N2D2L30) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        YAML::Node node = YAML::LoadFile("tests/w_1.0_D_2_N_2_L_30.yaml");
        node["test"] = true;
        node["jastrow"] = false;
        node["progress"] = false;
        node["importance"] = false;
        node["stepmc"] = 1.3;
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
       
        node["stepmc"] = 0.01;
        node["importance"] = true;
        E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30
    
    TEST(energiesw028N2D2L30) {
        /* check that energy in case without jastrow reproduces HF energy
         * (within std) */
        YAML::Node node = YAML::LoadFile("tests/w_0.28_D_2_N_2_L_20.yaml");
        node["test"] = true;
        node["jastrow"] = false;
        node["progress"] = false;
        node["importance"] = false;
        node["stepmc"] = 1.3;
        double E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
       
        node["stepmc"] = 2.5;
        node["importance"] = true;
        E = run(node);
        CHECK_CLOSE(node["E0"].as<double>(), E, 0.01);
    } // end TEST energiesw1N2D2L30
} // end SUITE VMCHF

#endif
