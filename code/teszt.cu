#include <iostream>
#include <math.h>
#include <chrono>
#include <iomanip>


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
void init(int len, float value, float* vector)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < len; i += stride)
        vector[i] = value;
}

__global__
void osszeg(int n, float *x, float *y)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
        for (int i = index; i < n; i += stride)
            y[i] += x[i];
}


__global__
void kiir(int n, int *a)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int osszegIndex = blockIdx.x;
    for(int i=index; i<n; i+=stride)
        a[i]=osszegIndex;
}

__global__
void toltessuruseg(int cellaSzam, int reszecskeSzam, int reszToltesSzam, int* indexek, float** reszRho)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int reszReszecskeSzam = (reszecskeSzam + reszToltesSzam - 1)/reszToltesSzam;
    for(int i=index; i<reszToltesSzam; i+=stride)
    {
        for(int j=i*reszReszecskeSzam; j<(i+1)*reszReszecskeSzam && j<reszecskeSzam; j++)
            reszRho[i][indexek[j]]+=1.0f;
    }
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

__global__
void besorol(int reszecskeSzam, int cellaSzam, float* helyek, int* indexek)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i=index; i<reszecskeSzam; i+=stride)
            indexek[i] = ((int)round(helyek[i])%cellaSzam + cellaSzam)%cellaSzam;
}

__global__ void ujE(int cellaSzam, float* fi1, float* fi2, float* eredmeny)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i=index; i<cellaSzam; i+=stride)
        eredmeny[i]=fi1[(i-1+cellaSzam)%cellaSzam]-fi2[(i+1)%cellaSzam];
}


int main(void)
{
    int Ng = 10, Nr = 11, Np = 200, i, j, *p;
    float **rRho, *rho, *sumRho;
    cudaMallocManaged(&rRho, Nr*sizeof(float*));
    for(i=0; i<Nr; i++)
        cudaMallocManaged(&rRho[i], Ng*sizeof(float));
    cudaMallocManaged(&p, Np*sizeof(int));
    cudaMallocManaged(&rho, Ng*sizeof(float));
    cudaMallocManaged(&sumRho, Ng*sizeof(float));

    init<<<4, 32>>>(Ng, 0, sumRho);
    cudaDeviceSynchronize();

    for(i=0; i<Np; i++)
    {
        p[i]=i%Ng;
        //cout << p[i] << endl;
    }

    //void toltessuruseg(int cellaSzam, int reszecskeSzam, int reszToltesSzam, int* indexek, float** reszRho)
    toltessuruseg<<<5, 32>>>(Ng, Np, Nr, p, rRho);
    cudaDeviceSynchronize();

    for(i=0; i<Nr; i++)
    {
        cout << "rRho[" << setw(2) << setfill(' ') << right << i << "]:  ";
        for(j=0; j<Ng; j++)
            cout << setw(4) << setfill(' ') << right << setprecision(4) << rRho[i][j];
        cout << endl;
    }

    for(i=0; i<Nr; i++)
        for(j=0; j<Ng; j++)
            sumRho[j]+=rRho[i][j];

    for(i=0; i<Np; i++)
        rho[p[i]]++;

    // a tényleg tesztelés alatt álló rész
    for(i=1; i<Nr; i++)
    {
        osszeg<<<5, 32>>>(Ng, rRho[i], rRho[0]);
        cudaDeviceSynchronize();
    }


    cout << endl << endl;

    cout << "rho =      ";
    for(j=0; j<Ng; j++)
        cout << setw(4) << setfill(' ') << right << setprecision(4) << rho[j];
    cout << endl;

    cout << "sumRho =   ";
    for(j=0; j<Ng; j++)
        cout << setw(4) << setfill(' ') << right << setprecision(4) << sumRho[j];
    cout << endl;

    cout << "rRho[0] =  ";
    for(j=0; j<Ng; j++)
        cout << setw(4) << setfill(' ') << right << setprecision(4) << rRho[0][j];
    cout << endl;


    // felszabadulás
    for(i=0; i<Nr; i++)
        cudaFree(rRho[i]);
    cudaFree(rRho);
    cudaFree(p);
    cudaFree(rho);
    cudaFree(sumRho);


    return 0;
}
