/*******************************************************************************************
 * RandomNumberGenerator.cpp
 * Functions to generate random numbers
 * Created by Gunnar Dreßler on 16.05.21.
 *******************************************************************************************/

#include "RandomNumberGenerator.h"

// declare the random number generator engine: there should only be a single engine
// used throughout the model!
mt19937 RandomNumberGenerator::rng_engine;

/*******************************************************************************************
 * Set the seed of the random number generator
 *******************************************************************************************/
void RandomNumberGenerator::set_random_seed(int seed) {
    // use the provided seed to initalize the engine
    get_random_engine().seed(seed);
    rng_engine.seed(seed);
}

void RandomNumberGenerator::set_random_seed() {
    // no random seed provided, use random_device to create a unique random seed
    random_device rd;
    get_random_engine().seed(rd());
}

/*******************************************************************************************
 * functions to generate random numbers
 * /todo: add generators for other distribtions, e.g. normal, exponential, etc.
 *******************************************************************************************/
double RandomNumberGenerator::generate_random_float(double min, double max) {
    uniform_real_distribution<double> uniform_dist(min, max);
    return uniform_dist(get_random_engine());
}

int RandomNumberGenerator::generate_random_int(int min, int max) {
    uniform_int_distribution<int> uniform_dist(min, max);
    return uniform_dist(get_random_engine());
}

double RandomNumberGenerator::generate_random_beta(double alpha, double beta){
    
    // if X is a random number drawn from gamma(alpha, 1) and Y is random number
    // drawn from gamma(beta, 1), then Z = X/(X+Y) is from beta(alpha, beta)

    gamma_distribution<double> dist_X(alpha, 1.0);
    gamma_distribution<double> dist_Y(beta, 1.0);

    double X = dist_X(get_random_engine());
    double Y = dist_Y(get_random_engine());

    return X / (X + Y);  
}

/*******************************************************************************************
 * return reference to the random number engine which is used in the generate_* functions
 *******************************************************************************************/
mt19937 &RandomNumberGenerator::get_random_engine() {

    return rng_engine;
}

/*******************************************************************************************
 * Constructor & destructor
 *******************************************************************************************/
RandomNumberGenerator::RandomNumberGenerator() = default;

RandomNumberGenerator::~RandomNumberGenerator() = default;
