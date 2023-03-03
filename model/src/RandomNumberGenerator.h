/*******************************************************************************************
 * RandomNumberGenerator.h
 * Functions to generate random numbers
 * Created by Gunnar Dreßler on 16.05.21.
 *******************************************************************************************/

#ifndef C___RANDOM_TEST_RANDOMNUMBERGENERATOR_H
#define C___RANDOM_TEST_RANDOMNUMBERGENERATOR_H

#include <random>

using namespace std;

class RandomNumberGenerator {

private:
    // Random number generator: Mersenne Twister engine to be used throughout the model
    static mt19937 rng_engine;
    // return reference to the random number engine
    static mt19937 &get_random_engine();

public:
    // constructor
    RandomNumberGenerator();

    // generate random numbers
    static double generate_random_float(double min, double max);
    static int generate_random_int(int min, int max);
    static double generate_random_beta(double alpha, double beta);
    
    // set random number generator seed
    static void set_random_seed(int seed);
    static void set_random_seed();

    // destructor
    virtual ~RandomNumberGenerator();
};

#endif //C___RANDOM_TEST_RANDOMNUMBERGENERATOR_H
