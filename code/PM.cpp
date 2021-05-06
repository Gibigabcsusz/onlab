#include <iostream>
#include <math.h>
#include <random>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <typeinfo>


using namespace std;
using namespace std::chrono;

void add(int n, float *x, float *y, float *z); //TODO
void ciklikus(float cellaSzam, int reszecskeSzam, float* helyek); //TODO
void besorol(int reszecskeSzam, int cellaSzam, float* helyek, int* indexek);
void filePrinter(int vektorHossz, float* x, float* y, string fileNev, string xLabel, string yLabel, string dataLabel);

int main(void)
{
    // Bemenetek megadása
    const int T = 1;
    const int Ta = 1; // az ábrázolás időlépésének száma, min=2
    const int Ng = 10000;
    const int Nc = 15;
    const int Np = Nc*Ng;
    const float maxvin = 0;
    const float omDT = 0.2;
    const float fihSzorzo = 100;
    const long int seedNum = 12;


    // Tároló vektorok inicializálása
    float *xa, *xb, **x, *va, *vb, **v, *rho, *fi, *fih, *Ea, *Eb, **E;
    int *p;

    xa = (float*)malloc(Np*sizeof(float));
    xb = (float*)malloc(Np*sizeof(float));
    x = (float**)malloc(2*sizeof(float*));
    va = (float*)malloc(Np*sizeof(float));
    vb = (float*)malloc(Np*sizeof(float));
    v = (float**)malloc(2*sizeof(float*));
    p = (int*)malloc(Np*sizeof(int));
    rho = (float*)malloc((Ng+1)*sizeof(float));
    fi = (float*)malloc((Ng+1)*sizeof(float));
    fih = (float*)malloc((Ng+1)*sizeof(float));
    Ea = (float*)malloc((Ng+1)*sizeof(float));
    Eb = (float*)malloc((Ng+1)*sizeof(float));
    E = (float**)malloc(2*sizeof(float*));

    x[0]=xa;
    x[1]=xb;
    v[0]=va;
    v[1]=vb;
    E[0]=Ea;
    E[1]=Eb;

    bool vOverkill = false;
    float vatlag = 0;
    int i, t;
    float rhoSzorzo = omDT*omDT/2/Nc;
    // TODO
    float y=0;
    for(i=0; i<Ng+1; i++)
    {
        fih[i] = fihSzorzo*(-256*y*y*y*y + 512*y*y*y - 352*y*y + 96*y - 9);
        y = y + 1.0/Ng;
    }
    float sorozat[Ng+1] = {};
    for(i=0; i<Ng+1; i++)
        sorozat[i]=i;

    // Kiinduló állapotok legenerálása
    uniform_real_distribution<float> unifx(-0.5, Ng-0.5);
    uniform_real_distribution<float> unifv(-maxvin, maxvin);
    default_random_engine re;
    re.seed(seedNum);
    for(i=0; i<Np; i++)
    {
        x[0][i] = unifx(re);
        v[0][i] = unifv(re);
        x[1][i] = x[0][i]+v[0][i];
    }
    ciklikus(Np, Ng, x[1]);

    // időmérés indítása
    auto check4 = high_resolution_clock::now();
    auto check5 = high_resolution_clock::now();
    auto check6 = high_resolution_clock::now();
    auto check7 = high_resolution_clock::now();
    auto check8 = high_resolution_clock::now();
    auto check9 = high_resolution_clock::now();
    auto startt = high_resolution_clock::now();

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
    //TODO
    for(i=0; i<Ng; i++)
    {
        fi[i]+=fih[i];
    }

    auto check2 = high_resolution_clock::now();

    // A térerősségek és egyben a gyorsulások ellentettjeinek kiszámolása
    E[0][0]=fi[Ng-1]-fi[1];
    E[0][Ng-1]=fi[Ng-2]-fi[0];
    //TODO
    for(i=1;i<Ng-1;i++)
    {
        E[0][i]=fi[i-1]-fi[i+1];
    }


    auto check3 = high_resolution_clock::now();


    // Átlagsebesség kiszámolása
    for(i=0;i<Np;i++)
        vatlag+=v[0][i];
    vatlag/=Np;




    ///////////////////
    // a nagy ciklus //
    ///////////////////

    for(t=1;t<T+1;t++)
    {

    check4 = high_resolution_clock::now();

        // a részecskéket cellákhoz rendeljük a legújabb kiszámolt pozíciók alapján
        besorol(Np, Ng, x[t%2], p);


    check5 = high_resolution_clock::now();

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
        // Ábra generálása
        /*if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng+1, sorozat, fih, "output/fih.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng+1, sorozat, fi, "output/fi.dat", "X", "Potenciál", "Pot");
            filePrinter(Ng+1, sorozat, rho, "output/rho.dat", "X", "Töltéssűrűség", "rho");
        }*/
        //TODO
        for(i=0; i<Ng; i++)
        {
            fi[i]+=fih[i]*fihSzorzo;
        }
        // Ábra generálása
        /*if(t==Ta)
        {
            fi[Ng]=fi[0];
            filePrinter(Ng+1, sorozat, fi, "output/fi2.dat", "X", "Potenciál", "Pot");
        }*/

    check6 = high_resolution_clock::now();

        // A térerősségek és egyben a gyorsulások ellentettjeinek kiszámolása
        E[t%2][0]=fi[Ng-1]-fi[1];
        E[t%2][Ng-1]=fi[Ng-2]-fi[0];
        //TODO
        for(i=1;i<Ng-1;i++)
        {
            E[t%2][i]=fi[i-1]-fi[i+1];
        }

    check7 = high_resolution_clock::now();

        // A következő sebességek kiszámolása
        for(i=0;i<Np;i++) //TODO
        // ez volt eredetileg:            v[t%2][i] = v[(t-1)%2][i] - (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
            v[t%2][i] = v[(t-1)%2][i] + (E[(t-1)%2][p[i]]+E[(t%2)][p[i]])/2;
        for(i=0;i<Np;i++)
        {
            vatlag += v[t%2][i];
            if(abs(v[t%2][i]) > Ng)
                vOverkill = true;
        }
        vatlag/=Np;
            // az átlagsebességgel korrigálás TODO
        for(i=0;i<Np;i++)
        {
            v[t%2][i]-=vatlag;
        }

    check8 = high_resolution_clock::now();

        // A következő helyek kiszámolása
        for(i=0;i<Np;i++)
            x[(t+1)%2][i] = x[t%2][i] + (v[(t-1)%2][i] + v[t%2][i])/2;

    check9 = high_resolution_clock::now();

        ciklikus(Ng, Np, x[(t+1)%2]);


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

    if(vOverkill)
        cout << "Túl nagy sebesség..." << endl;

    cout << "Sta-Sto: " << fullmicrosecs << " ms" << endl;
    cout << "Sta-Ch1: " << microsecs1 << " ms" << endl;
    cout << "Ch1-Ch2: " << microsecs2 << " ms" << endl;
    cout << "Ch2-Ch3: " << microsecs3 << " ms" << endl;
    cout << "Ch3-Ch4: " << microsecs4 << " ms" << endl;
    cout << "Ch4-Ch5: " << microsecs5 << " ms" << endl;
    cout << "Ch5-Ch6: " << microsecs6 << " ms" << endl;
    cout << "Ch6-Ch7: " << microsecs7 << " ms" << endl;
    cout << "Ch7-Ch8: " << microsecs8 << " ms" << endl;
    cout << "Ch8-Ch9: " << microsecs9 << " ms" << endl;
    cout << "Ch9-Sto: " << microsecs10 << " ms" << endl;


    //felszabadítás

    free(xa);
    free(xb);
    free(x);
    free(va);
    free(vb);
    free(v);
    free(p);
    free(rho);
    free(fi);
    free(fih);
    free(Ea);
    free(Eb);
    free(E);

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
void ciklikus(float cellaSzam, int reszecskeSzam, float* helyek) //TODO
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
    myFile.close();
}
