#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <typeinfo>


using namespace std;
using namespace std::chrono;

__global__ void osszegAXBY(int n, float a, float b, float *x, float *y, float *z);
__global__ void osszeg(int n, float *x, float *y);
__global__ void ujE(int cellaSzam, float* fi1, float* fi2, float* eredmeny);
__global__ void ciklikus(float cellaSzam, int reszecskeSzam, float* helyek);
__global__ void besorol(int reszecskeSzam, int cellaSzam, float* helyek, int* indexek);
__global__ void init(int len, float value, float* vector);
__global__ void toltessuruseg(int cellaSzam, int reszecskeSzam, int reszToltesSzam, int* indexek, float** reszRho);

void fihGenerator(int len, float szorzo, float* outVector);
void filePrinter(int vektorHossz, float* x, float* y, string fileNev, string xLabel, string yLabel, string dataLabel);
void kezdetiXV(long int seed, float maxv, int reszecskeSzam, int cellaSzam, float** xp, float** vp);

int main(void)
{
    // Bemenetek megadása
    const int T = 200;
    const int Ta = 199; // az ábrázolás időlépésének száma, min=2
    const int Ng = 1000;
    const int Nc = 15;
    const int Np = Nc*Ng;
    const float maxvin = 1;
    const float omDT = 0.2;
    const float fihSzorzo = 1;
    const long int seedNum = 0;
    const int blockSize = 32;
    const int numBlocksGrid = (Ng + blockSize - 1) / blockSize;
    const int numBlocksParticles = (Np + blockSize - 1) / blockSize;
    const int Nr = 1<<3; // résztöltések száma, érdemes 2^n alakúnak lennie


    float *sorozat = (float*)malloc(Ng*sizeof(float));
    float rhoSzorzo = omDT*omDT/2/Nc;
    int i, t;
    float vatlag = 0;

    // Tároló vektorok inicializálása a device-on
    float **x, **v, *rho, *fi, *fih, **E, *fiMasolat, **rRho, *egysegGrid;
    int *p;

    cudaMallocManaged(&x, 2*sizeof(float*));
    cudaMallocManaged(&egysegGrid, Ng*sizeof(float*));
    cudaMallocManaged(&v, 2*sizeof(float*));
    cudaMallocManaged(&(x[0]), Np*sizeof(float));
    cudaMallocManaged(&(x[1]), Np*sizeof(float));
    cudaMallocManaged(&(v[0]), Np*sizeof(float));
    cudaMallocManaged(&(v[1]), Np*sizeof(float));
    cudaMallocManaged(&p, Np*sizeof(int));
    cudaMallocManaged(&rho, Ng*sizeof(float));
    cudaMallocManaged(&fi, Ng*sizeof(float));
    cudaMallocManaged(&fiMasolat, Ng*sizeof(float));
    cudaMallocManaged(&fih, Ng*sizeof(float));
    cudaMallocManaged(&E, 2*sizeof(float*));
    cudaMallocManaged(&(E[0]), Np*sizeof(float));
    cudaMallocManaged(&(E[1]), Np*sizeof(float));
    cudaMallocManaged(&rRho, Nr*sizeof(float*));
    for(i=0; i<Nr; i++)
        cudaMallocManaged(&rRho[i], Ng*sizeof(float));

    // egységvektor generálás
    init<<<numBlocksGrid, blockSize>>>(Ng, 1, egysegGrid);
    cudaDeviceSynchronize();

    // a háttérpotenciál vektorának generálása
    fihGenerator(Ng, fihSzorzo, fih);

    // az ábrázoláshoz használt vektor generálása
    for(i=0; i<Ng; i++)
        sorozat[i]=i;

    // Kiinduló állapotok legenerálása
    kezdetiXV(seedNum, maxvin, Np, Ng, x, v);
    ciklikus<<<blockSize, numBlocksParticles>>>((float)Ng, Np, x[1]);
    cudaDeviceSynchronize();

    // időmérés indítása
    auto startt = high_resolution_clock::now();

    // Részecskék cellákba sorolása
    besorol<<<blockSize, numBlocksParticles>>>(Np, Ng, x[0], p);
    cudaDeviceSynchronize();

    // a részecskék besorolása, innen a rho. Ennek skálázása
    init<<<numBlocksGrid, blockSize>>>(Ng, -Nc, rho);
    cudaDeviceSynchronize();
    for(i=0; i<Nr; i++)
    {
        init<<<numBlocksGrid, blockSize>>>(Ng, 0, rRho[i]);
        cudaDeviceSynchronize();
    }
    toltessuruseg<<<numBlocksGrid, blockSize>>>(Ng, Np, Nr, p, rRho);
    cudaDeviceSynchronize();
    for(i=0; i<Nr; i++)
    {
        osszeg<<<numBlocksGrid, blockSize>>>(Ng, rRho[i], rho);
        cudaDeviceSynchronize();
    }
    osszegAXBY<<<numBlocksGrid, blockSize>>>(Ng, rhoSzorzo, 0.0f, rho, rho, rho);
    cudaDeviceSynchronize();

    // Potenciál kiszámolása
    fi[0]=0;
    for(i=0; i<Ng; i++)
        fi[0]+=(i+1)*rho[i];
    fi[0]/=Ng;
    fi[Ng-1]=0;
    fi[1]=rho[0]+2*fi[0];
    for(i=2; i<Ng-1; i++)
    {
        fi[i]=rho[i-1]+2*fi[i-1]-fi[i-2];
    }
    osszegAXBY<<<blockSize, numBlocksGrid>>>(Ng, fihSzorzo, 1, fih, fi, fi);
    cudaDeviceSynchronize();

    // A térerősségek és egyben a gyorsulások ellentettjeinek kiszámolása
    cudaMemcpy(fiMasolat, fi, Ng*sizeof(float), cudaMemcpyDeviceToDevice);
    ujE<<<blockSize, numBlocksGrid>>>(Ng, fi, fiMasolat, E[0]);
    cudaDeviceSynchronize();

    // Átlagsebesség kiszámolása TODO
    for(i=0;i<Np;i++)
        vatlag+=v[0][i];
    vatlag/=Np;


    ///////////////////
    // a nagy ciklus //
    ///////////////////

    for(t=1;t<T+1;t++)
    {
        // a részecskéket cellákhoz rendeljük a legújabb kiszámolt pozíciók alapján
        besorol<<<blockSize, numBlocksParticles>>>(Np, Ng, x[t%2], p);
        cudaDeviceSynchronize();


        // a járulékos töltéssűrűségének kiszámolása, innen a rho. Ennek skálázása
        init<<<numBlocksGrid, blockSize>>>(Ng, -Nc, rho);
        cudaDeviceSynchronize();
        for(i=0; i<Nr; i++)
        {
            init<<<numBlocksGrid, blockSize>>>(Ng, 0, rRho[i]);
            cudaDeviceSynchronize();
        }
        toltessuruseg<<<numBlocksGrid, blockSize>>>(Ng, Np, Nr, p, rRho);
        cudaDeviceSynchronize();
        for(i=0; i<Nr; i++)
        {
            osszeg<<<numBlocksGrid, blockSize>>>(Ng, rRho[i], rho);
            cudaDeviceSynchronize();
        }
        osszegAXBY<<<numBlocksGrid, blockSize>>>(Ng, rhoSzorzo, 0.0f, rho, rho, rho);
        cudaDeviceSynchronize();


        // Potenciál kiszámolása EZ NEM MEGY GYORSABBAN
        fi[0]=0;
        for(i=0; i<Ng; i++)
            fi[0]+=(i+1)*rho[i];
        fi[0]/=Ng;
        fi[Ng-1]=0;
        fi[1]=rho[0]+2*fi[0];
        for(i=2; i<Ng-1; i++)
        {
            fi[i]=rho[i-1]+2*fi[i-1]-fi[i-2];
        }
        for(i=0; i<Ng; i++)

        // Ábra generálása
        if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng, sorozat, fih, "output/fih.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng, sorozat, fi, "output/fi.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng, sorozat, rho, "output/rho.dat", "X", "Töltéssűrűség", "rho");
        }
        // külső potenciál hozzáadása
        osszegAXBY<<<blockSize, numBlocksGrid>>>(Ng, fihSzorzo, 1, fih, fi, fi);
        cudaDeviceSynchronize();

        // Ábra generálása...
        if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng, sorozat, fi, "output/fi2.dat", "X", "Potenciál", "Pot");
        }

        // A térerősség (ami a gyorsulás ellentettje) kiszámolása
        cudaMemcpy(fiMasolat, fi, Ng*sizeof(float), cudaMemcpyDeviceToDevice);
        ujE<<<blockSize, numBlocksGrid>>>(Ng, fi, fiMasolat, E[t%2]);
        cudaDeviceSynchronize();

        // A következő sebességek kiszámolása EZ NEM MEGY GYORSABBAN, HACSAK NEM LEHET EGYSZERRE UGYANONNAN OLVASNI
        for(i=0;i<Np;i++) //TODO
        // ez volt eredetileg:            v[t%2][i] = v[(t-1)%2][i] - (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
            v[t%2][i] = v[(t-1)%2][i] + (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
        for(i=0;i<Np;i++) //TODO
        {
            vatlag += v[t%2][i];
            if(abs(v[t%2][i]) > Ng)
                cout << "Túl nagy sebesség!" << endl;
        }
        vatlag/=Np;
            // az átlagsebességgel korrigálás TODO
        for(i=0;i<Np;i++)
        {
            v[t%2][i]-=vatlag;
        }



        // A következő helyek kiszámolása TODO
        for(i=0;i<Np;i++)
            x[(t+1)%2][i] = x[t%2][i] + (v[(t-1)%2][i] + v[t%2][i])/2;

        ciklikus<<<blockSize, numBlocksParticles>>>(Ng, Np, x[(t+1)%2]);
        cudaDeviceSynchronize();
    }

    auto stopp = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(stopp - startt);
    int millisecs = duration.count();
    cout << "Chrono meres: " << millisecs << " ms" << endl;




    // felszabadítás
    free(sorozat);

    cudaFree(egysegGrid);
    cudaFree(x[0]);
    cudaFree(x[1]);
    cudaFree(x);
    cudaFree(v[0]);
    cudaFree(v[1]);
    cudaFree(v);
    cudaFree(p);
    cudaFree(rho);
    cudaFree(fi);
    cudaFree(fih);
    cudaFree(E[0]);
    cudaFree(E[1]);
    cudaFree(E);
    for(i=0; i<Nr; i++)
        cudaFree(rRho[i]);
    cudaFree(rRho);




    return 0;
}


// function to osszegAXBY the elements of two arrays
__global__
void osszegAXBY(int n, float a, float b, float *x, float *y, float *z)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if(b==0)
        for (int i = index; i < n; i += stride)
            z[i] = a*x[i];
    else
        for (int i = index; i < n; i += stride)
            z[i] = a*x[i] + b*y[i];
}


// egy egyszerűbb és gyorsabb függvény két n hosszú vektor
// elemenkénti összeadására,
// aminek az eredménye a második vektorban tárolódik el
__global__ void osszeg(int n, float *x, float *y)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < n; i += stride)
        y[i] = x[i] + y[i];
}


// ez a ciklikus szimulációs teret biztosítja, egyelőre csak akkor, ha
// a részecskék legnagyobb sebessége kisebb, mint Ng
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



// a helyek változóban tárolt részecske x értékeket besorolja a cellákba és a cellák indexeit eltárolja
__global__
void besorol(int reszecskeSzam, int cellaSzam, float* helyek, int* indexek)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i=index; i<reszecskeSzam; i+=stride)
        indexek[i] = ((int)round(helyek[i])%cellaSzam + cellaSzam)%cellaSzam;
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
void ujE(int cellaSzam, float* fi1, float* fi2, float* eredmeny)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for(int i=index; i<cellaSzam; i+=stride)
        eredmeny[i]=fi1[(i-1+cellaSzam)%cellaSzam]-fi2[(i+1)%cellaSzam];
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



void fihGenerator(int len, float szorzo, float* outVector)
{
    float y=0;
    for(int i=0; i<len; i++)
    {
        outVector[i] = szorzo*(-256*y*y*y*y + 512*y*y*y - 352*y*y + 96*y - 9);
        y = y + 1.0/len;
    }
}

void filePrinter(int vektorHossz, float* x, float* y, string fileNev, string xLabel, string yLabel, string dataLabel)
{
    ofstream myFile;
    myFile.open(fileNev);
    myFile << "# " << xLabel << " " << yLabel << endl;
    for(int i=0; i<vektorHossz; i++)
        myFile << x[i] << " " << y[i] << endl;
    myFile << vektorHossz << " " << y[0] << endl;
    myFile.close();
}

void kezdetiXV(long int seed, float maxv, int reszecskeSzam, int cellaSzam, float** xp, float** vp)
{
    uniform_real_distribution<float> unifx(-0.5, cellaSzam-0.5);
    uniform_real_distribution<float> unifv(-maxv, maxv);
    default_random_engine re;
    re.seed(seed);
    for(int i=0; i<reszecskeSzam; i++)
    {
        xp[0][i] = unifx(re);
        vp[0][i] = unifv(re);
        xp[1][i] = xp[0][i]+vp[0][i];
    }
}
