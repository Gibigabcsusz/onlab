#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>


using namespace std;

// function to add the elements of two arrays TODO
void add(int n, double *x, double *y, double *z)
{
    for (int i = 0; i < n; i++)
        z[i] = x[i] + y[i];
}


void ciklikus(int cellaSzam)
{
    
}

// TODO a helyek változóban tárolt részecske x értékeket besorolja a cellákba és a cellák indexeit eltárolja
void besorol(int reszecskeSzam, int cellaSzam, double* helyek, int* indexek)
{
    for(int k=0; k<reszecskeSzam; k++)
            indexek[k] = ((int)round(helyek[k])%cellaSzam + cellaSzam)%cellaSzam;
}


int main(void)
{
    // Bemenetek megadása
    const int T = 3001;
    const int Ng = 20;
    const int Nc = 15;
    const int Np = Nc*Ng;
    const double maxvin = 5;
    const double omDT = 0.2;
    const double fihSzorzo = 0.2;
    const long int seedNum = 4;


    // Tároló vektorok inicializálása
    double xa[Np] = {};
    double xb[Np] = {};
    double* x[2] = {xa, xb};
    double va[Np] = {};
    double vb[Np] = {};
    double* v[2] = {va, vb};
    int p[Np] = {};
    double rho[Ng] = {};
    double fi[Ng] = {};
    double fih[Ng] = {};
    double E[Ng] = {};
    int i;
    double y=0;
    // TODO
    for(i=0; i<Ng; i++)
    {
        fih[i] = 256*y*y*y*y - 512*y*y*y + 352*y*y - 96*y + 9;
        y = y + 1.0/Ng;
        //cout << "fih[" << i << "]:  " << fih[i] << endl;
    }


    // Kiinduló állapotok legenerálása
    uniform_real_distribution<double> unif(-0.5, Ng-0.5);
    default_random_engine re;
    re.seed(seedNum);
    for(i=0; i<Np; i++)
    {
        x[0][i] = unif(re);
        x[1][i] = unif(re);
        v[0][i] = x[1][i]-x[0][i];
    }


    // Részecskék cellákba sorolása
    besorol(Np, Ng, x[0], p);


    // Töltéssűrűség kiszámolása
    for(i=0; i<Ng; i++)
        rho[i] = -Nc;



    return 0;
}
