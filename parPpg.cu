#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <stdint.h>
#include <inttypes.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

#include "cycleTimer.h"
#include "utils.h"

#define DEBUG           1 
#if DEBUG
#define cudaCheckError(ans)  cudaAssert((ans), __FILE__, __LINE__);
inline void cudaAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr, "CUDA Error: %s at %s:%d\n",
        cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#else
#define cudaCheckError(ans) ans
#endif

#define METHOD          2
#define USE_FIRST_PR    1
#define SKIP_EVEN       1

#define THREADS_PER_BLK 128
#define SHARED_MEM_SIZE 8192
#define NUM_STREAMS     16
#define NUM_FIRST_PR    512 


//void print_primes(bool *masks, uint64_t lo, 
//                uint64_t n, bool valForPrime)

#if METHOD == 1
__global__ void soe(bool *masks, uint64_t n, uint64_t k)
{
    // Ignore idx = 0 and 1 so add 2
    uint64_t idx = blockIdx.x * blockDim.x + threadIdx.x + 2;

    // Mark 0 and 1 as non primes
    if (idx == 2) {
        masks[0] = true;
        masks[1] = true;
    }

    if (idx <= k) {
        // If current number is unmarked
        if (!masks[idx]) {
            // Mark all multiples of i
            // Start from 2i since we don't want to mark i
            for(int j=2*idx; j<=n; j+=idx) {
                masks[j] = true;
            } 
        }
    }
}
#elif METHOD == 2
#if USE_FIRST_PR
__global__ void soe(bool *masks, bool *primeMasks, uint64_t lo,
        uint64_t chunkLo, uint64_t n, uint64_t k)
#else
__global__ void soe(bool *masks, uint64_t lo, uint64_t chunkLo,
        uint64_t n, uint64_t k)
#endif
{
    __shared__ bool cache[SHARED_MEM_SIZE];

#if USE_FIRST_PR
    __shared__ bool firstPrimes[NUM_FIRST_PR+1];
#endif

    uint64_t tid = threadIdx.x;
    uint64_t off = blockIdx.x * blockDim.x;

    uint64_t blockPos = (blockIdx.x * SHARED_MEM_SIZE);
    uint64_t startPos = chunkLo + blockPos;
    uint64_t endPos = startPos + SHARED_MEM_SIZE;
    endPos = (n < endPos) ? n+1 : endPos;

    //printf("%lu %" PRIu64 " %" PRIu64 "\n", blockIdx.x, startPos, endPos); 
    uint64_t i, j;

    // Load respective portion from array
    for (i=startPos+tid; i<endPos; i+=THREADS_PER_BLK) {
        cache[i-startPos] = masks[i - chunkLo];
#if USE_FIRST_PR
        // Init firstPrimes array upto SHARED_MEM_SIZE
        //firstPrimes[i-startPos] = false;
#endif
        if (i - chunkLo > CHUNK_SIZE) {
            printf("LD cl:%"PRIu64",sp:%"PRIu64"\n",chunkLo, startPos);
        }
    }
#if USE_FIRST_PR
    // Initialize rest of the indices for firstPrimes
    for (i=tid; i<=NUM_FIRST_PR; i+=THREADS_PER_BLK) {
        firstPrimes[i] = primeMasks[i];
    }
#endif
    __syncthreads();

    // Only first thread should mark 0 and 1 as non primes
    // This happens when startPos = 0 and off+tid = 0
    if (!(startPos | (off + tid))) {
        cache[0] = true;
#if SKIP_EVEN == 0
        cache[1] = true;
#endif
    }

#if SKIP_EVEN
    i = (2*tid) + 3;
    for ( ; i<=k; i+= 2*THREADS_PER_BLK) {
#else
    for (i=tid+2; i<=k; i+= THREADS_PER_BLK) {
#endif
#if USE_FIRST_PR
        // If current number is unmarked
        // Start from 2 since 0 and 1 are not primes
        // Chedk only till the size of the array
#if SKIP_EVEN
        if ((i > NUM_FIRST_PR) || (!firstPrimes[i/2])) {
#else
        if ((i > NUM_FIRST_PR) || (!firstPrimes[i])) {
#endif
#endif
#if SKIP_EVEN
            // Find first multiple of i within the range
            uint64_t firstMult = (((2*startPos) + 1) / i) * i;
            // Need to get first odd since indices for even don't exist
            firstMult = (firstMult % 2) ? firstMult : firstMult + i;
            // and convert it back to correct range
            firstMult /= 2;//(firstMult) ? ((firstMult-1)/2) : (i/2);
#else
            // Find first multiple of i within the range
            uint64_t firstMult = (startPos / i) * i;
#endif
            // Mark all multiples of i within the range
            for(j=firstMult; j<endPos; j+=i) {
                // Do not mark i itself
#if SKIP_EVEN
                // Compare converted value of i with j
                if ((i/2 != j) && (j >= startPos)) {
#else
                if ((i != j) && (j >= startPos)) {
#endif
                    cache[j-startPos] = true;
                    //printf("%lu, %lu, s: %lu, e:%lu\n", i, j, startPos, endPos);
#if 0
                    if (j <= NUM_FIRST_PR) {
                        firstPrimes[j-2] = true;
                    }
#endif
                    if (j-startPos > SHARED_MEM_SIZE) {
                        printf("%"PRIu64"\n",j-startPos);
                    }
                }
            }   
#if USE_FIRST_PR
        }   
#endif
    }

    __syncthreads();

    // Copy back to global memory
    for (i=startPos + tid; i<endPos; i+=THREADS_PER_BLK) {
        masks[i - chunkLo] = cache[i-startPos];
        if (i - chunkLo > CHUNK_SIZE) {
            printf("ST cl:%"PRIu64",sp:%"PRIu64"\n",chunkLo, startPos);
            //printf("ST %"PRIu64"\n",startPos - chunkLo + i);
        }
    }
}
#endif

__global__ void print_arr(bool *arr, uint64_t lo, uint64_t chunkLo, uint64_t hi
        , bool cmp) 
{
#if SKIP_EVEN
    if((chunkLo <= 2) && (hi >= 2)) {
        // Print 2
        printf("2\n");
    }
#endif
#if 0
    // Find first odd number close to chunkLo
    uint64_t j = (chunkLo % 2) ? chunkLo : chunkLo + 1;
    for (; j<=hi; j+=2) {
#else
    for (uint64_t j=chunkLo; j<=hi; j++) {
#endif
        if(arr[j-chunkLo] == cmp) {
#if SKIP_EVEN
            printf("%"PRIu64"\n", (2*j) + 1); 
#else
            printf("%"PRIu64"\n", j);
#endif
        }   
    }   
}

/*void init_dev_masks(bool *devMasks, int n)
{
    // We should be able to mark elements upto n
    // Therefore, index can go upto n, requiring n+1 elements
    size_t maskSize = sizeof(bool) * (n + 1); 
    cudaMalloc(&devMasks, maskSize);
    cudaMemset(devMasks, 0, maskSize);
}

void deinit_dev_masks(bool *devMasks)
{
    cudaFree(devMasks);
}*/
#if USE_FIRST_PR
__global__ void soa(bool *masks, bool *primeMasks, uint64_t lo, uint64_t chunkLo,
        uint64_t n, uint64_t k)
#else
__global__ void soa(bool *masks, uint64_t lo, uint64_t chunkLo, uint64_t n,
        uint64_t k)
#endif
{
    __shared__ int cache[SHARED_MEM_SIZE];
    
    // LUTs for calculating modulo values
    __shared__ bool arr1[MOD_VAL];
    __shared__ bool arr2[MOD_VAL];
    __shared__ bool arr3[MOD_VAL];
    
#if USE_FIRST_PR
    __shared__ bool firstPrimes[NUM_FIRST_PR+1];
#endif

    uint64_t tid = threadIdx.x;
    uint64_t off = blockIdx.x * blockDim.x;

    uint64_t blockPos = (blockIdx.x * SHARED_MEM_SIZE);
    uint64_t startPos = chunkLo + blockPos;
    uint64_t endPos = startPos + SHARED_MEM_SIZE;
    endPos = (n < endPos) ? n+1 : endPos;

    //printf("%lu %" PRIu64 " %" PRIu64 "\n", blockIdx.x, startPos, endPos); 
    uint64_t i;

    // Load respective portion from array
    for (i=startPos+tid; i<endPos; i+=THREADS_PER_BLK) {
        cache[i-startPos] = (masks[i - chunkLo]) ? 1 : 0;
        if (i - chunkLo > CHUNK_SIZE) {
            printf("LD cl:%"PRIu64",sp:%"PRIu64"\n",chunkLo, startPos);
        }   
    }   

#if USE_FIRST_PR
    // Initialize rest of the indices for firstPrimes
    for (i=tid; i<=NUM_FIRST_PR; i+=THREADS_PER_BLK) {
        firstPrimes[i] = primeMasks[i];
    }   
#endif

    for (i=tid; i<MOD_VAL; i+=THREADS_PER_BLK) {
        arr1[i] = false;
        arr2[i] = false;
        arr3[i] = false;
    }

    __syncthreads();

    if (tid == 0) {
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
    }
    
    if (!((tid + off) | startPos)) {  
        // Mark 2, 3 and 5 as primes
        cache[2] = 1;
        cache[3] = 1;
        cache[5] = 1;
    }
    __syncthreads();

    // For solutions to 4x^2 + y^2 = n within the range,
    // max x and y are as follows
    uint64_t xmax = ceilf(sqrtf(endPos)/2);
    uint64_t ymax = ceilf(sqrtf(endPos));
    uint64_t x = 1;
    uint64_t y = 1;
    uint64_t x2 = 1;
    uint64_t y2 = 1;
    uint64_t pdt;
    uint8_t r;

    // Solutions for 4x^2 + y^2 = n
    for (x=tid+1; x*x<endPos; x+=THREADS_PER_BLK) {
        x2 = x * x;
        for (y=1; y*y<endPos; y+=2) {
            y2 = y * y;
            pdt = (4*x2) + y2;
            r = pdt % 60;
            if (pdt >= endPos) {
                // Break if this happens as y is increasing so 
                // pdt will keep increasing for same x
                break;
            }
            if (pdt >= startPos) {
                if(arr1[r]) {
                    // Flip bit
                   //cache[pdt-startPos] ^= true;
                    atomicXor(&cache[pdt-startPos], 1);
                }
            }
        }
    }
    
    // Solutions for 3x^2 + y^2
    xmax = ceilf(sqrtf(endPos/3));
   
    // Odd xs even ys
    for (x=((2*tid)+1); x*x<endPos; x+=THREADS_PER_BLK) {
        x2 = x * x;
        for (y=2; y*y<endPos; y+=2) {
            y2 = y * y;
            pdt = (3*x2) + y2; 
            r = pdt % 60; 
            if (pdt >= endPos) {
                // Break if this happens as y is increasing so 
                // pdt will keep increasing for same x
                break;
            }   
            if (pdt >= startPos) {
                if(arr2[r]) {
                    // Flip bit
                    atomicXor(&cache[pdt-startPos], 1);
                }   
            }   
        }   
    }
    
    // Solutions for 3x^2 - y^2
    // Max x will be when y = x because x > y condition has to be met
    xmax = ceilf(sqrtf(endPos/2));
    
    for (x=tid+1; x*x<endPos; x+=THREADS_PER_BLK) {
        x2 = x * x;
        for (i=1; i<x; i+=2) {
            y = x - i; // Get odd/even combos
            y2 = y * y;
            pdt = (3*x2) - y2;
            r = pdt % 60;
            if (pdt >= endPos) {
                // Break if this happens as y is decreasing so 
                // pdt will keep increasing for same x
                break;
            }
            if (pdt >= startPos) {
                if(arr3[r]) {
                    // Flip bit
                    atomicXor(&cache[pdt-startPos], 1);
                    //cache[pdt-startPos] ^= true;
                }
            }
        }
    }
    
    __syncthreads();

    // Mark multiples of squares
    uint64_t i2 = 4;
    uint64_t imax = ceilf(sqrtf(endPos)); 

    for (i=tid+2; i*i<endPos; i+=THREADS_PER_BLK) {
        i2 = i * i;
        
#if USE_FIRST_PR
        if ((i > NUM_FIRST_PR) || (!firstPrimes[i])) {
#endif
            // Find first multiple within range
            uint64_t firstMult = (startPos / i2) * i2;

            // If prime, mark all multiples of square of that prime
            for (uint64_t j=firstMult; j < endPos; j += i2) {
                if (j >= startPos) {
                    cache[j-startPos] = 0;
                }
            }
#if USE_FIRST_PR
        }
#endif
    }

    __syncthreads();

    // Copy back to global memory
    for (i=startPos + tid; i<endPos; i+=THREADS_PER_BLK) {
        masks[i - chunkLo] = (cache[i-startPos]) ? true : false;
        if (i - chunkLo > CHUNK_SIZE) {
            printf("ST cl:%"PRIu64",sp:%"PRIu64"\n",chunkLo, startPos);
            //printf("ST %"PRIu64"\n",startPos - chunkLo + i);
        }
    }

}

void par_algo(bool *devMasks, bool *primeMasks, uint64_t lo, uint64_t chunkLo,
        uint64_t n, uint64_t k, cudaStream_t stream, bool isSOA)
{
    //uint64_t k = floor(sqrt(n));
#if METHOD == 1
    int numBlks = (k / THREADS_PER_BLK) + 1;
#elif METHOD == 2
    int numBlks = (CHUNK_SIZE / SHARED_MEM_SIZE) + 1;
    if (n - chunkLo < CHUNK_SIZE) {
        numBlks = ((n - chunkLo) / SHARED_MEM_SIZE) + 1;
    }
#endif
    if (isSOA) {
#if USE_FIRST_PR
        soa<<<numBlks, THREADS_PER_BLK, 0, stream>>>(devMasks, primeMasks, lo,
                chunkLo, n, k);
#else
        soa<<<numBlks, THREADS_PER_BLK, 0, stream>>>(devMasks, lo, chunkLo, n, k);
#endif
    } else {
#if METHOD == 1
    soe<<<numBlks, THREADS_PER_BLK, 0, stream>>>(devMasks, n, k);
#else
#if USE_FIRST_PR
        soe<<<numBlks, THREADS_PER_BLK, 0, stream>>>(devMasks, primeMasks, lo,
                chunkLo, n, k);
#else
        soe<<<numBlks, THREADS_PER_BLK, 0, stream>>>(devMasks, lo, chunkLo, n, k);
#endif
#endif
    }
}

void par_find_primes(bool *hostMasks, uint64_t lo, uint64_t n, uint64_t k, bool isSOA,
        bool isDisplayTime, bool isQuiet, bool primeBoolVal)
{
    bool *devMasks[NUM_STREAMS];
    cudaStream_t streams[NUM_STREAMS]; 
    size_t maskSize = sizeof(bool) * (CHUNK_SIZE + 1); 
  
    bool *primeMasks;

    // Effective n (for skipping evens will be different)
#if SKIP_EVEN
    uint64_t nEff = (n-1)/2;
#else
    uint64_t nEff = n;
#endif

#if USE_FIRST_PR
    cudaMalloc(&primeMasks, sizeof(bool) * (NUM_FIRST_PR+1));
    cudaMemset(primeMasks, 0, sizeof(bool) * (NUM_FIRST_PR+1));
#endif

    // We should be able to mark elements upto n
    // Therefore, index can go upto n, requiring n+1 elements
    for (int i=0; i<NUM_STREAMS; i++) {
        cudaMalloc(&devMasks[i], maskSize);
        cudaMemset(devMasks[i], 0, maskSize);
        cudaCheckError(cudaStreamCreate(&streams[i]));
    }
    // Performance metrics for parallel execution
    cudaEvent_t start, stop;//, stop_mem;
    float duration;//, duration_mem;

    if (isDisplayTime) {
        // Init objects
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        //cudaEventCreate(&stop_mem);

        // Start recording time
        cudaEventRecord(start);
    }

#if METHOD == 2
#if USE_FIRST_PR
    soe<<<(NUM_FIRST_PR/CHUNK_SIZE) + 1, THREADS_PER_BLK>>>(primeMasks, 
            primeMasks, 0, 0, NUM_FIRST_PR, ceilf(sqrtf(NUM_FIRST_PR)));
    //print_arr<<<1,1>>>(primeMasks, 0, 0, NUM_FIRST_PR, false);
    //cudaDeviceSynchronize();
#endif
#endif
    int ctr = 0;
    for (uint64_t i=lo; i<=nEff; i+=CHUNK_SIZE) {
        uint64_t endPos = i + CHUNK_SIZE;
        endPos = (nEff < endPos) ? nEff : endPos;
        
        int streamId = ctr % NUM_STREAMS;
        cudaStream_t stream = streams[streamId];
        //printf("Stream %d\n", streamId);
        par_algo(devMasks[streamId], primeMasks, lo, i, endPos, k, stream, isSOA);
        
        if (!isQuiet) {
            print_arr<<<1, 1, 0, stream>>>(devMasks[streamId], lo, i, endPos, primeBoolVal);
        }
        ctr++;
        if (ctr >= NUM_STREAMS) {
            cudaCheckError(cudaStreamSynchronize(streams[ctr % NUM_STREAMS]));
            cudaCheckError(cudaMemsetAsync(devMasks[ctr % NUM_STREAMS], 0,
                        maskSize, streams[ctr % NUM_STREAMS]));
        }
    }
    
    if (isDisplayTime) {
        cudaEventRecord(stop);
    }

    //cudaMemcpy(hostMasks, devMasks, maskSize, cudaMemcpyDeviceToHost);
    
    if (isDisplayTime) {
        //cudaEventRecord(stop_mem);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&duration, start, stop);
        printf("Time taken by GPU for generating primes upto %lu: %fms\n",
                n, duration);

        //cudaEventSynchronize(stop_mem);
        //cudaEventElapsedTime(&duration_mem, start, stop_mem);
        //printf("Time taken by GPU for generating primes upto %lu, including"
        //        " memcopy: %fms\n", n, duration_mem);
    }

    for (int i=0; i<NUM_STREAMS; i++) {
        cudaStreamDestroy(streams[i]);
        cudaFree(devMasks[i]);
    }
}
