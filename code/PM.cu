#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>
#include <fstream>
#include <chrono>
//#include <typeinfo>


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
    const float fihSzorzo = 10;
    const long int seedNum = 12;
    const int blockSize = 32;
    const int numBlocksGrid = (Ng + blockSize - 1) / blockSize;
    const int numBlocksParticles = (Np + blockSize - 1) / blockSize;
    const int Nr = 8; // résztöltések száma, érdemes 2^n alakúnak lennie
    float *sorozat = (float*)malloc(Ng*sizeof(float));
    float rhoSzorzo = omDT*omDT/2/Nc;
    int i, t;
    float vatlag = 0;

    /////////////////
    // előkészítés //
    /////////////////

    // Tároló vektorok inicializálása a device-on
    float **x, **v, *rho, *fi, *fih, **E, *fiMasolat, **rRho;
    int *p;

    cudaMallocManaged(&x, 2*sizeof(float*));
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
    cudaMallocManaged(&(E[0]), Ng*sizeof(float));
    cudaMallocManaged(&(E[1]), Ng*sizeof(float));
    cudaMallocManaged(&rRho, Nr*sizeof(float*));
    for(i=0; i<Nr; i++)
        cudaMallocManaged(&rRho[i], Ng*sizeof(float));

    // a külső potenciál vektorának generálása és skálázása
    fihGenerator(Ng, fihSzorzo, fih);

    // az ábrázoláshoz használt vektor generálása, aminek az i-edik eleme i
    for(i=0; i<Ng; i++)
        sorozat[i]=i;

    // Kiinduló részecskesebességek és -helyek generálása
    kezdetiXV(seedNum, maxvin, Np, Ng, x, v);
    ciklikus<<<blockSize, numBlocksParticles>>>((float)Ng, Np, x[1]);
    cudaDeviceSynchronize();

    // időmérés indítása és a későbbi checkpointok deklarálása
    auto check4 = high_resolution_clock::now();
    auto check5 = high_resolution_clock::now();
    auto check6 = high_resolution_clock::now();
    auto check7 = high_resolution_clock::now();
    auto check8 = high_resolution_clock::now();
    auto check9 = high_resolution_clock::now();
    auto startt = high_resolution_clock::now();

    // Részecskék helyének cellákhoz rendelése
    besorol<<<blockSize, numBlocksParticles>>>(Np, Ng, x[0], p);
    cudaDeviceSynchronize();

    // A cellákhoz rendelt részecskék töltéssűrűségeinek rész-töltéssűrűségekre
    // bontása és ezek összegzése, hogy megkapjuk az össz-töltéssűrűséget.
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

auto check1 = high_resolution_clock::now();

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

auto check2 = high_resolution_clock::now();

    // A térerősségek és egyben a gyorsulások ellentettjeinek kiszámolása
    cudaMemcpy(fiMasolat, fi, Ng*sizeof(float), cudaMemcpyDeviceToDevice);
    ujE<<<blockSize, numBlocksGrid>>>(Ng, fi, fiMasolat, E[0]);
    cudaDeviceSynchronize();

auto check3 = high_resolution_clock::now();

    ///////////////////
    // a nagy ciklus //
    ///////////////////
    for(t=1;t<T+1;t++)
    {
    check4 = high_resolution_clock::now();

        // a részecskéket cellákhoz rendeljük a legújabb kiszámolt pozíciók alapján
        besorol<<<blockSize, numBlocksParticles>>>(Np, Ng, x[t%2], p);
        cudaDeviceSynchronize();

    check5 = high_resolution_clock::now();

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

    check6 = high_resolution_clock::now();

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

        // Ábra generálása
        if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng, sorozat, fih, "output/CUfih.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng, sorozat, fi, "output/CUfi.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng, sorozat, rho, "output/CUrho.dat", "X", "Töltéssűrűség", "rho");
        }
        // külső potenciál hozzáadása
        osszegAXBY<<<blockSize, numBlocksGrid>>>(Ng, fihSzorzo, 1, fih, fi, fi);
        cudaDeviceSynchronize();



        // Ábra generálása...
        if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng, sorozat, fi, "output/CUfi2.dat", "X", "Potenciál", "Pot");
        }

    check7 = high_resolution_clock::now();


        // A térerősség (ami a gyorsulás ellentettje) kiszámolása
        cudaMemcpy(fiMasolat, fi, Ng*sizeof(float), cudaMemcpyDeviceToDevice);
        ujE<<<blockSize, numBlocksGrid>>>(Ng, fi, fiMasolat, E[t%2]);
        cudaDeviceSynchronize();

    check8 = high_resolution_clock::now();

        // A következő sebességek kiszámolása EZ NEM MEGY GYORSABBAN, HACSAK NEM LEHET EGYSZERRE UGYANONNAN OLVASNI
        for(i=0;i<Np;i++) //TODO
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

    check9 = high_resolution_clock::now();

        // A következő helyek kiszámolása TODO
        for(i=0;i<Np;i++)
            x[(t+1)%2][i] = x[t%2][i] + (v[(t-1)%2][i] + v[t%2][i])/2;
        ciklikus<<<blockSize, numBlocksParticles>>>(Ng, Np, x[(t+1)%2]);
        cudaDeviceSynchronize();
    }

    auto stopp = high_resolution_clock::now();

    auto fullduration = duration_cast<microseconds>(stopp - startt);
    int fullmicrosecs = fullduration.count();
    auto duration1 = duration_cast<microseconds>(check1 - startt);
    int microsecs1 = duration1.count();
    auto duration2 = duration_cast<microseconds>(check2 - check1);
    int microsecs2 = duration2.count();
    auto duration3 = duration_cast<microseconds>(check3 - check2);
    int microsecs3 = duration3.count();
    auto duration4 = duration_cast<microseconds>(check4 - check3);
    int microsecs4 = duration4.count();
    auto duration5 = duration_cast<microseconds>(check5 - check4);
    int microsecs5 = duration5.count();
    auto duration6 = duration_cast<microseconds>(check6 - check5);
    int microsecs6 = duration6.count();
    auto duration7 = duration_cast<microseconds>(check7 - check6);
    int microsecs7 = duration7.count();
    auto duration8 = duration_cast<microseconds>(check8 - check7);
    int microsecs8 = duration8.count();
    auto duration9 = duration_cast<microseconds>(check9 - check8);
    int microsecs9 = duration9.count();
    auto duration10 = duration_cast<microseconds>(stopp - check9);
    int microsecs10 = duration10.count();
    cout << "Sta-Sto: " << fullmicrosecs << " us" << endl;
    cout << "Sta-Ch1: " << microsecs1 << " us" << endl;
    cout << "Ch1-Ch2: " << microsecs2 << " us" << endl;
    cout << "Ch2-Ch3: " << microsecs3 << " us" << endl;
    cout << "Ch3-Ch4: " << microsecs4 << " us" << endl;
    cout << "Ch4-Ch5: " << microsecs5 << " us" << endl;
    cout << "Ch5-Ch6: " << microsecs6 << " us" << endl;
    cout << "Ch6-Ch7: " << microsecs7 << " us" << endl;
    cout << "Ch7-Ch8: " << microsecs8 << " us" << endl;
    cout << "Ch8-Ch9: " << microsecs9 << " us" << endl;
    cout << "Ch9-Sto: " << microsecs10 << " us" << endl;

    // felszabadítás
    free(sorozat);
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
