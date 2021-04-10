#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>
#include <fstream>


using namespace std;

// function to add the elements of two arrays TODO
void add(int n, double *x, double *y, double *z)
{
    for (int i = 0; i < n; i++)
        z[i] = x[i] + y[i];
}


// ez a ciklikus szimulációs teret biztosítja, egyelőre csak akkor, ha
// a részecskék legnagyobb sebessége kisebb, mint Ng
void ciklikus(double cellaSzam, int reszecskeSzam, double* helyek) //TODO
{
    for(int k=0; k<reszecskeSzam; k++)
    {
        if(helyek[k] < -0.5)
            helyek[k] += cellaSzam;
        if(helyek[k] > cellaSzam-0.5)
            helyek[k] -= cellaSzam;
    }
}

// TODO a helyek változóban tárolt részecske x értékeket besorolja a cellákba és a cellák indexeit eltárolja
void besorol(int reszecskeSzam, int cellaSzam, double* helyek, int* indexek)
{
    for(int k=0; k<reszecskeSzam; k++)
            indexek[k] = ((int)round(helyek[k])%cellaSzam + cellaSzam)%cellaSzam;
}


void filePrinter(int vektorHossz, double* x, double* y, string fileNev)
{
    ofstream myFile;
    myFile.open(fileNev);
    for(int i=0; i<vektorHossz; i++)
        myFile << x[i] << " " << y[i] << endl;
    myFile.close();
}


int main(void)
{
    // Bemenetek megadása
    const int T = 3;
    const int Ng = 200;
    const int Nc = 15;
    const int Np = Nc*Ng;
    const double maxvin = 1;
    const double omDT = 0.2;
    const double fihSzorzo = 0.1;
    const long int seedNum = 6;


    // Tároló vektorok inicializálása
    double xa[Np] = {};
    double xb[Np] = {};
    double* x[2] = {xa, xb};
    double va[Np] = {};
    double vb[Np] = {};
    double* v[2] = {va, vb};
    int p[Np] = {};
    double rho[Ng+1] = {};
    double fi[Ng+1] = {};
    double fih[Ng+1] = {};
    double Ea[Ng+1] = {};
    double Eb[Ng+1] = {};
    double* E[2] = {Ea, Eb};
    double vatlag = 0;
    int i, t;
    double rhoSzorzo = omDT*omDT/2/Nc;
    // TODO
    double y=0;
    for(i=0; i<Ng; i++)
    {
        fih[i] = 256*y*y*y*y - 512*y*y*y + 352*y*y - 96*y + 9;
        y = y + 1.0/Ng;
    }


    // Kiinduló állapotok legenerálása
    uniform_real_distribution<double> unifx(-0.5, Ng-0.5);
    uniform_real_distribution<double> unifv(-maxvin, maxvin);
    default_random_engine re;
    re.seed(seedNum);
    for(i=0; i<Np; i++)
    {

        x[0][i] = unifx(re);
        v[0][i] = unifv(re);
        x[1][i] = x[0][i]+v[0][i];
    }
    ciklikus(Np, Ng, x[1]);


    // Részecskék cellákba sorolása
    besorol(Np, Ng, x[0], p);


    // Töltéssűrűség kiszámolása TODO
    for(i=0; i<Ng; i++) // A háttér-töltésűrűségre inicializálás
        rho[i] = -Nc;
    for(i=0; i<Np; i++) // A részecskék járulékos töltéssűrűségeinek hozzáadása
        rho[p[i]]++;
    int n = 9;
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
        fi[i]+=fih[i]*fihSzorzo;
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

    for(t=1;t<T;t++)
    {
        // a részecskéket cellákhoz rendeljük a legújabb kiszámolt pozíciók alapján
        besorol(Np, Ng, x[t%2], p);

        // Töltéssűrűség kiszámolása TODO
        for(i=0; i<Ng; i++) // A háttér-töltésűrűségre inicializálás
            rho[i] = -Nc;
        for(i=0; i<Np; i++) // A részecskék járulékos töltéssűrűségeinek hozzáadása
            rho[p[i]]++;
        int n = 9;
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
        //TODO
        for(i=0; i<Ng; i++)
        {
            fi[i]+=fih[i]*fihSzorzo;
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
            v[t%2][i] = v[(t-1)%2][i] - (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
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

//        cout << "Sebességek, t=" << t << endl;
//        for(i=0;i<Np;i++)
//        {
//            cout << fixed << setw(n) << setfill(' ') << setprecision(n-3) << right << v[t%2][i] << endl;
//        }


        // A következő helyek kiszámolása
        for(i=0;i<Np;i++)
            x[(t+1)%2][i] = x[t%2][i] + (v[(t-1)%2][i] + v[t%2][i])/2;
        ciklikus(Ng, Np, x[(t+1)%2]);


    }

    double sorozat[Ng+1] = {};
    for(i=0; i<Ng+1; i++)
        sorozat[i]=i;
    fi[Ng]=fi[0];
    filePrinter(Ng+1, sorozat, fi, "output/semmi.dat");


    return 0;
}
