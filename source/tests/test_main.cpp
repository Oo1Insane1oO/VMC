#ifdef HARMONICOSCILLATOR 
    #include "../tests/test_vmcqd.cpp"
    #include "../tests/test_jastrowho.cpp"
#endif

#ifdef HARTREEFOCK 
    #include "../tests/test_vmchf.cpp"
#endif

#include "../tests/test_slater.cpp"
#include "../tests/test_cartesian.cpp"
#include "../tests/test_methods.cpp"
#include "../tests/test_MTLS.cpp"


int test_main() {
    return UnitTest::RunAllTests();
} // end function test_main
