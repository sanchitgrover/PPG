#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#ifndef _UTILS_H_
#define _UTILS_H_

#define CHUNK_SIZE  10000000 
#define MOD_VAL     60

void print_primes(bool *masks, uint64_t lo, uint64_t n, bool valForPrime);

#endif
