#include <iostream> // to prevent length macro conflict
#include <Rinternals.h>

extern "C" {
void trait_model(double *input,double *stat_to_return);
}

