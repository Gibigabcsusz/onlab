#include <iostream>
#include <math.h>
#include <chrono>


using namespace std;
using namespace std::chrono;

// CUDA kernel to add elements of two arrays
__global__
void add(int n, float a, float b, float *x, float *y, float *z)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
        for (int i = index; i < n; i += stride)
            z[i] = a*x[i] + b*y[i];
}

__global__
void ciklikus(float cellaSzam, int reszecskeSzam, float* helyek)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i = index; i < reszecskeSzam; i += stride)
    {
        if(helyek[i] < -0.5)
            helyek[i] += cellaSzam;
        if(helyek[i] > cellaSzam-0.5)
            helyek[i] -= cellaSzam;
    }
}

int main(void)
{
    int N = 1<<5;
    float *x, *y;


    // Allocate Unified Memory -- accessible from CPU or GPU
    cudaMallocManaged(&x, N*sizeof(float));
    cudaMallocManaged(&y, N*sizeof(float));

    // initialize x and y arrays on the host
    for (int i = 0; i < N; i++) {
        x[i] = 1.0f;
        y[i] = 5.0f*i-N;
    }

    // Launch kernel on 1M elements on the GPU
    //int blockSize = 32;

    int blockSize = 32;
    int numBlocks = (N + blockSize - 1) / blockSize;

    auto startt = high_resolution_clock::now();

    //add<<<numBlocks, blockSize>>>(N, 1.0f, 1.0f, x, y, y);

    ciklikus<<<numBlocks, blockSize >>>(N, N, y);

    // Wait for GPU to finish before accessing on host
    cudaDeviceSynchronize();

    auto stopp = high_resolution_clock::now();

    auto duration = duration_cast<nanoseconds>(stopp - startt);
    int nanosecs = duration.count()%1000;
    int microsecs = ((duration.count()-nanosecs)/1000)%1000;
    int millisecs = ((duration.count()-nanosecs)/1000-microsecs)/1000;
    cout << "Chrono meres: " << millisecs << " ms  +  " << microsecs << " us  +  " << nanosecs << " ns" << endl;

    // Check for errors (all values should be 3.0f)
    float maxError = 0.0f;
    for (int i = 0; i < N; i++)
        cout << i << " - " << y[i] << endl;
    //    maxError = fmax(maxError, fabs(y[i]-3.0f));
    std::cout << "Max error: " << maxError << std::endl;


    // Free memory
    cudaFree(x);
    cudaFree(y);

    return 0;
}
