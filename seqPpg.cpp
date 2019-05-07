#include <cmath>
#include <vector>
#include <getopt.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "cycleTimer.h"
#include "utils.h"

#define SEQ_METHOD  2

void par_find_primes(bool *masks, uint64_t lo, uint64_t n, uint64_t k, bool isSOA, 
        bool isDisplayTime, bool isQuiet, bool primeVal);

void usage(const char* progname) {
    printf("Usage: %s [options]\n", progname);
    printf("Program Options:\n");
    printf("  -a  --algo <0/1>           Algorithm to be used. 0 for Sieve of Eratosthenes, 1 for Sieve of Atkin\n");
    printf("  -c  --check                Check correctness of output\n");
    printf("  -d  --display              Display the execution time of algorithm\n");
    printf("  -l  --lower  <LOWER>       Display prime numbers starting from this positive integer\n");
    printf("  -n  --limit  <LIMIT>       A positive integer to display prime numbers upto\n");
    printf("  -p  --parallel             Type of execution. 0 for sequential, 1 for parallel\n");
    printf("  -q  --quiet                If used, do not isplay the primes generated after execution\n");
    printf("  -?  --help                 This message\n");
}

void seq_soe(bool *masks, uint64_t n, bool *primeMasks,
        uint64_t lo, uint64_t k)
{
    if (lo == 0) {
        // Mark 0 and 1 as non-prime
        primeMasks[0] = true;
        masks[0] = true;
        primeMasks[1] = true;
        masks[1] = true;
    }
    
    //uint64_t endPos = (n < lo + CHUNK_SIZE) ? n : lo + CHUNK_SIZE;

    for (uint64_t i=2; i<=k; i++) {
        // Do following if current number is unmarked
        if (!primeMasks[i]) {
            // First multiple within the range
            uint64_t firstMult = (lo / i) * i;
            // Mark all multiples of i
            // Start from 2i since we don't want to mark i
            for(uint64_t j=firstMult; j<=n; j+=i) {
                if ((i != j) && (j >= lo)) {
                    masks[j-lo] = true;
                    // Mark in primeMasks also if j <= k
                    if (j <= k) {
                        primeMasks[j] = true;
                    }
                }
            }
        }
    }
}

void print_primes(bool *masks, uint64_t lo, uint64_t n, bool valForPrime)
{
    for (uint64_t i=lo; i<=n; i++) {
        if (masks[i-lo] == valForPrime) {
            printf("%lu\n", i);
        }
    }
}

void seq_soa(bool *masks, uint64_t n, bool *primeMasks, 
        uint64_t lo, uint64_t k)
{
    // LUTs for calculating modulo values
    bool arr1[MOD_VAL] = {false};
    bool arr2[MOD_VAL] = {false};
    bool arr3[MOD_VAL] = {false};
    
    // Set appropriate indices to true
    arr1[1] = true;
    arr1[13] = true;
    arr1[17] = true;
    arr1[29] = true;
    arr1[37] = true;
    arr1[41] = true;
    arr1[49] = true;
    arr1[53] = true;

    arr2[7] = true;
    arr2[19] = true;
    arr2[31] = true;
    arr2[43] = true;

    arr3[11] = true;
    arr3[23] = true;
    arr3[47] = true;
    arr3[59] = true;

    // Mark 2, 3 and 5 as primes
    masks[2] = true;
    masks[3] = true;
    masks[5] = true;

    // Initial square of x and y for value 1
    uint64_t x = 1;
    uint64_t y = 1;
    uint64_t x2 = 1;
    uint64_t y2 = 1;

#if SEQ_METHOD == 1
    for (; x2 < n; x++) {
        x2 = x * x;
        y2 = 1;
        for (y=1; y2 < n; y++) {
            y2 = y * y;

            unsigned int pdt1 = (4 * x2) + y2;
            unsigned int pdt2 = (3 * x2) + y2;
            unsigned int pdt3 = (3 * x2) - y2;
            // Checking for mod60 remainders
            unsigned int r1 = pdt1 % MOD_VAL;
            unsigned int r2 = pdt2 % MOD_VAL;
            unsigned int r3 = pdt3 % MOD_VAL;

            // Mod60 remainder for this product should be 1, 13, 17, 
            // 29, 37, 41, 49 or 53
            if ((pdt1 <= n) && (arr1[r1])) {
                // Flip entry for pdt
                masks[pdt1] = !masks[pdt1];
            }

            // Mod60 remainder for this product should be 7, 19, 31 or 43
            if ((pdt2 <= n) && (arr2[r2])) {
                 masks[pdt2] = !masks[pdt2];
            }

            // Mod60 remainder for this product should be 11, 23, 47 or 59
            if ((x > y) && (pdt3 <= n) && (arr3[r3])) {
                 masks[pdt3] = !masks[pdt3];
            }
        }
    }
#elif SEQ_METHOD == 2
    for (; x2 < n; x++) {
        x2 = x * x;
        y2 = 1;
        for (y=1; y2 < n; y+=2) {
            y2 = y * y;

            unsigned int pdt = (4 * x2) + y2; 
            unsigned int r = pdt % MOD_VAL;
            if ((pdt <= n) && (arr1[r])) {
                masks[pdt] = !masks[pdt];
            }
        }
    }
    
    x2 = 1;
    for (x=1; x2 < n; x+=2) {
        x2 = x * x;
        y2 = 1;
        for (y=2; y2 < n; y+=2) {
            y2 = y * y;
            unsigned int pdt = (3 * x2) + y2;
            unsigned int r = pdt % MOD_VAL;
        
            if ((pdt <= n) && (arr2[r])) {
                masks[pdt] = !masks[pdt];
            }
        }
    }

    x2 = 1;
    for (x=2; x2 < n; x++) {
        x2 = x * x;
        for (uint64_t i=1; i < x; i+=2) {
            y = x - i;
            y2 = y * y;
            unsigned int pdt = (3 * x2) - y2;
            unsigned int r = pdt % MOD_VAL;
        
            if ((pdt <= n) && (arr3[r])) {
                masks[pdt] = !masks[pdt];
            }
        }
    }
#endif
    // Precompute first i^2 for use later
    uint64_t i2 = 4; // For i = 2
    
    // Find the next prime number in list
    for (uint64_t i=2; i2 < n; i++) {
        i2 = i * i;
        if ((i>=lo) && (masks[i])) {
            // If prime, mark all multiples of square of that prime
            for (uint64_t j=i2; j < n; j += i2) {
                if (j >= lo) {
                    masks[j] = false;
                }
            }
        }
    }
}

void init_masks(bool *masks, uint64_t n)
{
    memset(masks, 0, (n+1) * sizeof(bool));
}

void seq_find_primes(bool *masks, uint64_t n, bool* primeMasks, 
        uint64_t k, bool isSOA, bool isDisplayTime, bool isQuiet, bool primeBoolVal)
{
    // Use SOE by default
    void (*func)(bool *, uint64_t, bool *, uint64_t,
            uint64_t);
    func  = &seq_soe;

    if (isSOA) {
        func = &seq_soa;
    }

    double startTime;
    if (isDisplayTime) {
        startTime = CycleTimer::currentSeconds();
    }
    for (uint64_t i=0; i<n; i+=CHUNK_SIZE) {
        uint64_t endPos = i + CHUNK_SIZE;
        endPos = (n < endPos) ? n : endPos;
        func(masks, endPos, primeMasks, i, k);
        if (!isQuiet) {
            print_primes(masks, i, endPos, primeBoolVal);
        }
        init_masks(masks, CHUNK_SIZE);
    }
    if (isDisplayTime) {
        double endTime = CycleTimer::currentSeconds() - startTime;
        printf("Time taken by CPU for generating primes upto %lu: %fs\n", n, endTime);
    }
}

int main(int argc, char **argv)
{
    bool checkCorrectness = false;

    // Number to find primes upto
    uint64_t n = 1000;

    // Lower bound on starting prime numbers
    uint64_t lo = 0;

    // Run sequential by default
    bool isParallel = false;

    // Run SOE by default
    bool isSOA = false;

    // Do not display execution time
    bool isDisplayTime = false;

    // Display prime numbers generated by default
    bool isQuiet = false;

    // Prime numbers are labelled false in masks array
    bool primeBoolVal = false;

    // parse commandline options ////////////////////////////////////////////
    int opt;
    static struct option long_options[] = {
        {"help",     0, 0,  '?'},
        {"check",    0, 0,  'c'},
        {"display",  0, 0,  'd'},
        {"lower",    1, 0,  'l'},
        {"limit",    1, 0,  'n'},
        {"parallel", 0, 0,  'p'},
        {"quiet", 0, 0,  'q'},
        {"algo",     1, 0,  'a'},
        {0 ,0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "a:l:n:pqcd?", long_options, NULL)) != EOF) {

        switch (opt) {
        case 'c':
            checkCorrectness = true;
            break;
        case 'd':
            isDisplayTime = true;
            break;
        case 'l':
            lo = strtoul(optarg, NULL, 10);
            break;
        case 'n':
            n = strtoul(optarg, NULL, 10);
            break;
        case 'p':
            isParallel = true;
            break;
        case 'q':
            isQuiet = true;
            break;
        case 'a':
            // If algo is 1, switch to SOA
            if (atoi(optarg) == 1) {
                isSOA = true;
                primeBoolVal = true;
            }
            break;
        case '?':
        default:
            usage(argv[0]);
            return 1;
        }
    }
    // end parsing of commandline options //////////////////////////////////////
    
    // Bit vector for masking numbers
    bool *masks = (bool *) malloc(sizeof(bool) * (CHUNK_SIZE + 1));

    // Init masks
    init_masks(masks, CHUNK_SIZE);

    // Store prime numbers upto sqrt(n)
    uint64_t k = floor(sqrt(n));

    bool *primeMasks = (bool *) malloc(sizeof(bool) * (k + 1));
    init_masks(primeMasks, k);

    if (isParallel) {
        par_find_primes(masks, lo, n, k, isSOA, isDisplayTime, isQuiet, primeBoolVal);
    } else {
        seq_find_primes(masks, n, primeMasks, k, isSOA, isDisplayTime, isQuiet,
                primeBoolVal);
    }

    if(checkCorrectness) {
        printf("Checking\n");
    }
    
    free(primeMasks);
    free(masks);
    return 0;
}
