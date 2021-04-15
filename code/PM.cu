#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>
#include <fstream>

/*    __global__
     void add(int n, float *x, float *y)
     {
       int index = blockIdx.x * blockDim.x + threadIdx.x;
       int stride = blockDim.x * gridDim.x;
       for (int i = index; i < n; i += stride)
         y[i] = x[i] + y[i];
     }
*/

using namespace std;

void add(int n, float a, float b, float *x, float *y, float *z); //TODO
void ciklikus(float cellaSzam, int reszecskeSzam, float* helyek); //TODO
void besorol(int reszecskeSzam, int cellaSzam, float* helyek, int* indexek);//TODO
void filePrinter(int vektorHossz, float* x, float* y, string fileNev, string xLabel, string yLabel, string dataLabel);
void fihGenerator(int len, float szorzo, float* outVector);//TODO
void kezdetiXV(long int seed, float maxv, int reszecskeSzam, int cellaSzam, float** xp, float** vp);
__global__ void init(int len, float value, float* vector);

int main(void)
{
    // Bemenetek megadása
    const int T = 200;
    const int Ta = 199; // az ábrázolás időlépésének száma, min=2
    const int Ng = 100;
    const int Nc = 15;
    const int Np = Nc*Ng;
    const float maxvin = 1;
    const float omDT = 0.2;
    const float fihSzorzo = 1;
    const long int seedNum = 0;


    // Tároló vektorok inicializálása a host-on
/*    float *x[2] = {(float*)malloc(Np*sizeof(float)), (float*)malloc(Np*sizeof(float))};
    float *v[2] = {(float*)malloc(Np*sizeof(float)), (float*)malloc(Np*sizeof(float))};
    int *p = (int*)malloc(Np*sizeof(int));
    float *rho = (float*)malloc(Ng*sizeof(float));
    float *fi = (float*)malloc(Ng*sizeof(float));
    float *fih = (float*)malloc(Ng*sizeof(float));
    float *E[2] = {(float*)malloc(Ng*sizeof(float)), (float*)malloc(Ng*sizeof(float))}; */
    float *sorozat = (float*)malloc(Ng*sizeof(float));
    float rhoSzorzo = omDT*omDT/2/Nc;
    int i, t;
    float vatlag = 0;

    // Tároló vektorok inicializálása a device-on
    // DUPLA HOSSZÚ VEKTOROK KELLENEK
    float **x, **v, *rho, *fi, *fih, **E;
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
    cudaMallocManaged(&fih, Ng*sizeof(float));
    cudaMallocManaged(&E, 2*sizeof(float*));
    cudaMallocManaged(&(E[0]), Np*sizeof(float));
    cudaMallocManaged(&(E[1]), Np*sizeof(float));




    fihGenerator(Ng, fihSzorzo, fih);

    for(i=0; i<Ng; i++)
        sorozat[i]=i;

    // Kiinduló állapotok legenerálása
    kezdetiXV(seedNum, maxvin, Np, Ng, x, v);

    // Részecskék cellákba sorolása
    besorol(Np, Ng, x[0], p);


    // Töltéssűrűség kiszámolása TODO
    for(i=0; i<Ng; i++) // A háttér-töltésűrűségre inicializálás
        rho[i] = -Nc;
    for(i=0; i<Np; i++) // A részecskék járulékos töltéssűrűségeinek hozzáadása
        rho[p[i]]++;
    for(i=0; i<Ng; i++) // Az eredmény skálázása
    {
        rho[i] *= rhoSzorzo;
    }

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
    //TODO
    for(i=0; i<Ng; i++)
    {
        fi[i]+=fih[i];
    }

    // A térerősségek és egyben a gyorsulások ellentettjeinek kiszámolása
    E[0][0]=fi[Ng-1]-fi[1];
    E[0][Ng-1]=fi[Ng-2]-fi[0];
    //TODO
    for(i=1;i<Ng-1;i++)
    {
        E[0][i]=fi[i-1]-fi[i+1];
    }

    // Átlagsebesség kiszámolása
    for(i=0;i<Np;i++)
        vatlag+=v[0][i];
    vatlag/=Np;




    ///////////////////
    // a nagy ciklus //
    ///////////////////

    for(t=1;t<T+1;t++)
    {
        // a részecskéket cellákhoz rendeljük a legújabb kiszámolt pozíciók alapján
        besorol(Np, Ng, x[t%2], p);

        // Töltéssűrűség kiszámolása TODO
        for(i=0; i<Ng; i++) // A háttér-töltésűrűségre inicializálás
            rho[i] = -Nc;
        for(i=0; i<Np; i++) // A részecskék járulékos töltéssűrűségeinek hozzáadása
            rho[p[i]]++;
        for(i=0; i<Ng; i++) // Az eredmény skálázása
        {
            rho[i] *= rhoSzorzo;
        }

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
        for(i=0; i<Ng; i++)
        // Ábra generálása
        if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng, sorozat, fih, "output/fih.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng, sorozat, fi, "output/fi.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng, sorozat, rho, "output/rho.dat", "X", "Töltéssűrűség", "rho");
        }
        //TODO
        for(i=0; i<Ng; i++)
        {
            fi[i]+=fih[i]*fihSzorzo;
        }
        // Ábra generálása
        if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng, sorozat, fi, "output/fi2.dat", "X", "Potenciál", "Pot");
        }
        // A térerősségek és egyben a gyorsulások ellentettjeinek kiszámolása
        E[t%2][0]=fi[Ng-1]-fi[1];
        E[t%2][Ng-1]=fi[Ng-2]-fi[0];
        //TODO
        for(i=1;i<Ng-1;i++)
        {
            E[t%2][i]=fi[i-1]-fi[i+1];
        }

        // A következő sebességek kiszámolása
        for(i=0;i<Np;i++) //TODO
        // ez volt eredetileg:            v[t%2][i] = v[(t-1)%2][i] - (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
            v[t%2][i] = v[(t-1)%2][i] + (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
        for(i=0;i<Np;i++)
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



        // A következő helyek kiszámolása
        for(i=0;i<Np;i++)
            x[(t+1)%2][i] = x[t%2][i] + (v[(t-1)%2][i] + v[t%2][i])/2;
        ciklikus(Ng, Np, x[(t+1)%2]);



    }


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




    return 0;
}


// function to add the elements of two arrays TODO
void add(int n, float *x, float *y, float *z)
{
    for (int i = 0; i < n; i++)
        z[i] = x[i] + y[i];
}


// ez a ciklikus szimulációs teret biztosítja, egyelőre csak akkor, ha
// a részecskék legnagyobb sebessége kisebb, mint Ng
void ciklikus(float cellaSzam, int reszecskeSzam, float* helyek)
{
    for(int k=0; k<reszecskeSzam; k++)
    {
        if(helyek[k] < -0.5)
            helyek[k] += cellaSzam;
        if(helyek[k] > cellaSzam-0.5)
            helyek[k] -= cellaSzam;
    }
}

// a helyek változóban tárolt részecske x értékeket besorolja a cellákba és a cellák indexeit eltárolja
void besorol(int reszecskeSzam, int cellaSzam, float* helyek, int* indexek)
{
    for(int k=0; k<reszecskeSzam; k++)
            indexek[k] = ((int)round(helyek[k])%cellaSzam + cellaSzam)%cellaSzam;
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

void fihGenerator(int len, float szorzo, float* outVector)
{
    float y=0;
    for(int i=0; i<len; i++)
    {
        outVector[i] = szorzo*(-256*y*y*y*y + 512*y*y*y - 352*y*y + 96*y - 9);
        y = y + 1.0/len;
    }
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
    ciklikus(reszecskeSzam, cellaSzam, xp[1]);
}

__global__
void init(int len, float value, float* vector)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    for (int i = index; i < len; i += stride)
        vector[i] = value;
}

