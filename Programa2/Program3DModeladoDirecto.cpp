
/* pi */
#include "pch.h"
#include <iostream>
#include<complex>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdlib.h>
#include <omp.h>
#include <iomanip>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl.h>

#define M_PI  3.14159265358979323846
struct IncGenerator {
	int current_;
	IncGenerator(int start) : current_(start) {}
	int operator() () { return current_++; }
};

using namespace std;
typedef struct {
	double re;
	double i;
}
doublecomplex;
vector<int>vecID, vecID3;
bool myfunction(int i, int j) { return (vecID[i] < vecID[j]); };
float Lx = 30.0, Ly = 30.0, Lz = 30.0; //Dimensiones del dominio en Km
int eleX, eleY, eleZ, num1z, num2z, num3z, num4z, num5z,num6z,q;
float *x, *y, *z, dx, dy, dz, z1, z2, x1, x2, dp;
int num1x, num2x, num3x, num4x, num5x;
int num1y, num2y, num3y, num4y, num5y;
double facX, Sum, facZ, facZ1, facZ2;

int nx, ny, nz,totS;
void CuerpoCircular();
void CuerpoHomogeneo();
void CuerpoRectangular();
void CuerpoCircular() {
	/* Definición del numero de elementos y nodos en el plano XY */
	num1x = 15;
	num2x = 16;
	num1y = 15;
	num2y = 16;
	eleX = (int)(2 * num1x + num2x);
	eleY = (int)(2 * num1x + num2x);
	num1z = 25;
	num2z = 10;
	num3z = 16;
	num4z = 20;
	eleZ = (int)(num1z + num2z + num3z + num4z);
	nx = eleX + 1;
	ny = eleY + 1;
	nz = eleZ + 1;	
	int p1,p2;
	p1 = 0;
	facX = 0.2; Sum = 0.2; facZ = 0.05; facZ1 = 1.2; facZ2 = 1.2;
	y = (float*)malloc(ny * sizeof(float));
	x = (float*)malloc(nx * sizeof(float));

	/* Definición de coordenadas en la dirección X */
	for (int i = 0; i < (num1x); i++)
	{
		p2 = num1x - i - 1;
		x[p2] = -(float)1.0 - (float)Sum;//-(float)(dx*(float)(totS/2.0)+(float)Sum);//
		p1++;
		facX = facX * 1.3;
		Sum = Sum + facX;
	}
	dx = 2.0 / (float)(num2x);
	for (int i = 0; i < (num2x + 1); i++)
	{
		x[p1] = -1.0 + (float)dx*i;
		p1++;
	}
	facX = 0.2;
	for (int i = 0; i < (num1x); i++)
	{
		x[p1] = x[p1 - 1] + (float)facX;
		p1++;
		facX = facX * 1.3;
	}

	/* Definición de coordenadas en la dirección Y */
	for (int i = 0; i < ny; i++)
	{
		y[i] = x[i];
	}
}
void CuerpoRectangular() {
	/* Definición del numero de elementos y nodos en el plano XY */
	num1x = 13; num2x = 6; num3x = 8;//num2x=6
	num1y = 13; num2y = 6; num3y = 8;
	eleX = (int)(2 * num1x + 2 * num2x + num3x);
	eleY = (int)(2 * num1x + 2 * num2x + num3x);
	num1z = 20;//38
	num2z = 10;
	num3z = 6;
	num4z = 8;
	num5z = 6;
	num6z = 15;//num5z=14;
	eleZ = (int)(num1z + num2z + num3z + num4z + num5z + num6z);
	cout << "eleX: " << eleX << "eleY: " << eleY << "eleZ: " << eleZ << endl;
	nx = eleX + 1;
	ny = eleY + 1;
	nz = eleZ + 1;	

	y = (float*)malloc(ny * sizeof(float));
	x = (float*)malloc(nx * sizeof(float));
	
	/*Definición de coordenadas en X*/
	int p1, p2;
	p1 = 0;
	facX = 0.2; Sum = 0.2; facZ = 0.05; facZ1 = 1.2; facZ2 = 1.2;
	dx = 1.5 / (float)num3x;
	totS = (int)(2 * num2x + num3x);
	for (int i = 0; i < (num1x); i++)
	{
		p2 = num1x - i - 1;
		x[p2] = -(float)1.0 - (float)Sum;//-(float)(dx*(float)(totS/2.0)+(float)Sum);//
		p1++;
		facX = facX * 1.3;
		Sum = Sum + facX;
	}
	dx = 0.5 / (double)num2x;
	for (int i = 0; i < (num2x + 1); i++)
	{
		x[p1] = -1.0 + (float)dx*i;
		p1++;
	}
	dx = 1.0 / (double)num3x;
	for (int i = 1; i < (num3x + 1); i++)
	{
		x[p1] = -0.5 + (float)dx*i;
		p1++;
	}
	dx = 0.5 / (double)num2x;
	for (int i = 1; i < (num2x + 1); i++)
	{
		x[p1] = 0.5 + (float)dx*i;
		p1++;
	}
	facX = 0.2;
	//Sum=0.2;
	for (int i = 0; i < (num1x); i++)
	{
		x[p1] = x[p1 - 1] + (float)facX;
		p1++;
		facX = facX * 1.3;
	}
	/*Definición de coordenadas en Y*/
	for (int i = 0; i < ny; i++)
	{
		y[i] = x[i];
	}
}
void CuerpoHomogeneo()
{
	/* Definición del numero de elementos y nodos en el plano XY */
	eleX = 45;//36;
	eleY = 45;//36;
	num1z = 25;//25;
	num2z = 38;//35
	eleZ = num1z + num2z;
	cout << "eleX: " << eleX << "eleY: " << eleY << "eleZ: " << eleZ << endl;
	nx = eleX + 1;
	ny = eleY + 1;
	nz = eleZ + 1;
	
	dx = Lx / (float)(nx - 1);
	dy = Ly / (float)(ny - 1);

	y = (float*)malloc(ny * sizeof(float));
	x = (float*)malloc(nx * sizeof(float));
	/*Definición de coordenadas en X*/
	for (int i = 0; i < nx; i++)
	{
		x[i] = -Lx / 2 + dx * (float)i;
	}
	/*Definición de coordenadas en Y*/
	for (int i = 0; i < ny; i++)
	{
		y[i] = -Ly / 2 + dy * (float)i;
		//cout << y[i] << endl;
	}
	int p1, p2;
	p1 = 0;
	facX = 0.2; Sum = 0.2; facZ = 0.05; facZ1 = 1.2; facZ2 = 1.2;

}
int CASO;
int main()
{	
	/* PARTE 1: Creación de malla rectilínea en el plano XY para cualquiera de 3 casos:
	   - CASO=1: Prisma rectangular embebido en medio homogéneo
	   - CASO=2: Medio homogéneo
	   - CASO=3: Cuerpo "esférico" embebido en medio homogéneo */
	CuerpoCircular(); //Caso de esfera
	CASO = 3;
	
	int skX, skY,skZ,ipos,jpos;	
	int aristasFronNum,p=0, elemtot, vartot;
	complex<float> H0(1.0, 0.0), A, B, C;
	double ResA;
	complex <float>c(2.0, 0.0);
	complex<float>val2(1000.0, 0.0);
	ResA = (double)3.3*pow(10.0, 8.0);	
	complex<double>k0, k1, val1;
	complex <double> im(0,1);
	double Freq = 100.0, wr, ResH = 100.0;
	complex<float>val3(0.0, 0.0);
	complex<double>valC(1000.0, 0.0);
	double mu0;
	mu0 = 4.0 * M_PI*pow(10.0, -7.0);
	C = H0;
	A.real((double)(H0.real() / 2.0 + H0.real()*pow(10.0, 4.0)*sqrt(3.3 / ResH) / 2.0));
	A.imag(0.0);
	B.real(C.real() - A.real());
	B.imag(C.imag() - A.imag());
	double valY;
	complex<float>varA(-1.0, 0.0);
	q = 0;	
	int q2;
	int chk;	
	int l[12];
	double Ke[12];
	double K1[4][4] = { {2.0,-2.0,1.0,-1.0},
					   {-2.0,2.0,-1.0,1.0},
					   {1.0,-1.0,2.0,-2.0},
					   {-1.0,1.0,-2.0,2.0} };
	double K2[4][4] = { {2.0,1.0,-2.0,-1.0},
					   {1.0,2.0,-1.0,-2.0},
					   {-2.0,-1.0,2.0,1.0},
					   {-1.0,-2.0,1.0,2.0} };
	double K3[4][4] = { {2.0,1.0,-2.0,-1.0},
					   {-2.0,-1.0,2.0,1.0},
					   {1.0,2.0,-1.0,-2.0},
					   {-1.0,-2.0,1.0,2.0} };
	double F[4][4] = { {4.0,2.0,2.0,1.0},
					   {2.0,4.0,1.0,2.0},
					   {2.0,1.0,4.0,2.0},
					   {1.0,2.0,2.0,4.0} };

	complex<double>K[12][12];
	double valP, valR, valI;
	int    numV;
	int valF, sizeF;
	int enume;
	complex<double>Sd;
	int id = 0, id2 = 0, max, ac = 0;
	int numS = 0;
	p2 = -1;
	doublecomplex refV;
	int numT;
	double a[2], b[2], A1[2], B1[2], A2[2], B2[2], C2[2];
	double RhoXY, RhoYX, RhoXX, RhoYY;
	double FaseYX, FaseXY;
	complex <double> Ezy, Ezx, Eyz, Exz, Zxx, Zxy, Zyx, Zyy, detZ;
	complex <double> *Ex1, *Ey1, *Hx1, *Hy1, *Ex2, *Ey2, *Hx2, *Hy2;
	int numT2 = (int)((eleX - 1)*(eleY - 1));
	Ex1 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Ey1 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Hx1 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Hy1 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Ex2 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Ey2 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Hx2 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);
	Hy2 = (complex<double>*)malloc(sizeof(complex<double>)*numT2);	
	
	/* Definición del arreglo de frecuencias a evaluar */
	double Fq[25];
	double f2, f1,df;
	f2 = log10(1000.0);
	f1 = log10(0.001);
	df = (f2 - f1) / (24.0);
	for (int i = 0; i < 25; i++)
	{
		Fq[i] = pow(10.0,f1+(double)(df*i));
	}
	//Fq[0] = 10;
	double limZ = Lz / 2.0;
	double ddz = log10(2.05)- log10(0.05);
	double dT,stZ=0;
	int numS2;
	double ResC=0.5,Sumac=0.0;
	double xmed, ymed, zmed, rdis;
	/*FOR principal que itera el valor de la frecuencia*/
	for (int fq = 0; fq < 1; fq++)
	{
		int *dirichI, *cond;
		double*z;
		vector <double> zP1, zP2, zP3;
		complex <float> *dirichV;
		int **numb;
		double*lx, *ly, *lz;
		double *Resis;
		complex<double> *res;
		int*ia;
		int*varc;
		int *IA, *JA;
		complex<double>*DV;
		doublecomplex*xR, *bv;
		Freq = Fq[fq];
		wr = 2.0 * M_PI*Freq;
		complex<double>vecR1(0.0, mu0 / ResA * wr), vecR2(0.0, mu0 / ResH * wr);
		k0 = sqrt(vecR1);
		k1 = sqrt(vecR2);
		/*PARTE 2: Redefinición del arreglo Z según el skin depth*/
		/*La distancia al límite inferior se define como 10 veces el skin depth*/
		limZ = 10.0*sqrt(2.0 * ResH / (mu0*wr)) / 1000.0;
		facZ = 0.001;
		Sumac = 0.0;
		double limI, limS;
		if (CASO == 3)
		{
			limI = 0.05;
			limS = 2.05;
		}
		if (CASO == 1)
		{
			limI = 0.1;
			limS = 2.0;
		}
		while (limZ > (Sumac + facZ))
		{
			Sumac = Sumac + facZ;
			if ((CASO == 1) || (CASO == 3))
			{
				if ((limI) > Sumac)//-num3z*dx
				{
					zP2.push_back(Sumac);
				}
				if ((limS) < Sumac)//+(double)num5z*dx
				{
					zP3.push_back(Sumac);
				}
			}
			zP1.push_back(Sumac);
			//cout<<Sumac<<endl;
			facZ = facZ * facZ1;
		}
		num1z = zP1.size() + 1;
		num2z = zP2.size() + 1;
		if (CASO == 1)
		{
			num6z = zP3.size() + 1;
			eleZ = num1z + num2z + num3z + num4z + num5z + num6z;
		}
		if (CASO == 3)
		{
			num4z = zP3.size() + 1;
			eleZ = num1z + num2z + num3z + num4z;
		}
		nz = eleZ + 1;
		z = (double*)malloc(nz * sizeof(double));
		cout << "num1z: " << num1z << endl;
		cout << "num2z: " << num2z << endl;
		//cout <<"num4z: "<<num4z << endl;
		cout << "eleZ: " << eleZ << endl;
		if (limZ <= (limS))//1.8+(double)num5z*dx
		{
			if (CASO == 3)
			{
				num4z = 0;
			}
			if (CASO == 1)
			{
				num6z = 0;
			}
		}

		z[0] = (double)limZ;
		p1 = 1;
		for (int i = 1; i < (num1z); i++)
		{
			p2 = num1z - 1 - i;
			z[p1] = zP1[p2];
			p1 = p1 + 1;
		}
		z[p1] = 0.0;
		p1++;
		for (int i = 1; i < (num2z); i++)
		{
			p2 = i - 1;
			z[p1] = -zP2[p2];
			//cout<<z[p1]<<endl;
			p1 = p1 + 1;
		}
		if (CASO == 3)
		{
			dx = 2.0 / (double)num3z;
		}
		if (CASO == 1)
		{
			dx = 0.4 / (double)num3z;
		}
		for (int i = 0; i < (num3z + 1); i++)
		{
			z[p1] = -infI - (double)dx*i;
			//cout<<z[p1]<<endl;
			p1 = p1 + 1;
		}
		if (CASO == 1)
		{
			dx = 1.1 / (double)num4z;
			for (int i = 1; i < (num4z); i++)
			{
				z[p1] = -0.5 - (double)dx*i;
				cout << z[p1] << endl;
				p1 = p1 + 1;
			}
			dx = 0.4 / (double)num5z;
			for (int i = 0; i < (num5z + 1); i++)
			{
				z[p1] = -1.6 - (double)dx*i;
				cout << z[p1] << endl;
				p1 = p1 + 1;
			}


			if (num6z != 0)
			{
				for (int i = 1; i < (num6z); i++)
				{
					p2 = i - 1;
					z[p1] = -zP3[p2];
					//cout<<z[p1]<<endl;
					p1 = p1 + 1;
				}
				z[p1] = -(double)limZ;
			}
		}
	

		if (CASO == 3)
		{
			if (num4z != 0)
			{
				for (int i = 1; i < (num4z); i++)
				{
					p2 = i - 1;
					z[p1] = -zP3[p2];
					//cout<<z[p1]<<endl;
					p1 = p1 + 1;
				}
				z[p1] = -(double)limZ;
			}
		}
		
		
		
		cout <<"Frecuencia: "<<Freq<<" Hz"<<endl;
		cout<<"limZ: " << limZ << endl;

		aristasFronNum = (int)((eleY * nz + eleZ * ny) * 2 + ((nx - 2)*eleZ + nz * eleX) * 2 + (eleY*(nx - 2) + eleX * (ny - 2)) * 2); //Número de aristas en la frontera
		elemtot = (int)(eleX*eleY*eleZ); //Número de prismas en todo el dominio
		vartot = (int)(((eleY*nz) + (eleZ)*ny)*nx + eleX * ny*nz);//Número de aristas en todo el dominio
		//dirichI = (int*)malloc(aristasFronNum * sizeof(int));
		
		cond = (int*)malloc(vartot * sizeof(int));	//Arreglo indicador de fronteras
	
		Resis = (double*)malloc(elemtot * sizeof(double));//Arreglo con valores de resistividad de cada celda	
	
		dirichV = (complex<float>*)malloc(vartot * sizeof(complex<float>));	 //Arreglo donde se guardaron los valores en las fronteras
		
		/*PARTE 3: Definición de los valores en las fronteras*/
		for (int i = 0; i < vartot; i++)
		{
			cond[i] = 0;
		}
				//Abajo dir X
		for (int j = 0; j < ny; j++)
		{
		for (int i = 0; i < eleX; i++)
		{
			p = (int)(eleY*nx*nz + eleZ * ny*nx) + (int)(i + eleX * j + (nz - 1)*eleX*ny);
			dirichV[p] = 0.0;
			cond[p] = 1;

		}
		}
		//Abajo dir Y

		for (int i = 0; i < nx; i++)
		{
		for (int j = 0; j < eleY; j++)
		{
			p = (int)(nz - 1)*(eleY*nx) + (int)(i*eleY) + j;
			dirichV[p] = 0.0;
			cond[p] = 1;

			q = q + 1;
		}
		}
	
		//Izquierda y derecha dir Z
		for (int i = 0; i < 2; i++)
		{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < eleZ; k++)
			{
				p = (int)(eleY*nx*nz) + (int)(k + eleZ * j + eleZ * ny*(nx - 1)*i);
				dirichV[p] = 0.0;
				//dirichI[q] = p + 1;
				q = q + 1;
				cond[p] = 1;
			}
		}
		}
		//Atrás y adelante dir Z
		for (int j = 0; j < 2; j++)
		{
		for (int i = 1; i < (nx - 1); i++)
		{
			for (int k = 0; k < eleZ; k++)
			{
				p = (int)(eleY*nx*nz) + (int)(k + eleZ * (ny - 1)*j + eleZ * ny*i);
				dirichV[p] = 0.0;
				//dirichI[q] = p + 1;
				q = q + 1;
				cond[p] = 1;
			}
		}
		}
		chk = (int)(eleY*nx*nz + eleZ * ny*nx + 10);
		cout << "VariablesNum:" << vartot << endl;
		cout << "FronterasNum:" << aristasFronNum << endl;
		lx = (double*)malloc(elemtot * sizeof(double));
		ly = (double*)malloc(elemtot * sizeof(double));
		lz = (double*)malloc(elemtot * sizeof(double));

		numb = (int**)malloc(elemtot * sizeof(int *));
	

		for (int i = 0; i < elemtot; i++)
		{
		numb[i] = (int*)malloc(12 * sizeof(int));
		}
		
		IA = (int*)malloc((vartot + 1) * sizeof(int));
		JA = (int*)malloc(vartot * sizeof(int));	
		DV = (complex<double>*)malloc(vartot * sizeof(complex<double>));		
		varc = (int*)malloc(vartot * sizeof(int));	
		numV = (int)(vartot - aristasFronNum);
		cout <<"Size of matrix:"<<numV << endl;
		xR = (doublecomplex*)malloc(numV * sizeof(doublecomplex));
		bv = (doublecomplex*)malloc(numV * sizeof(doublecomplex));
		ia = (int*)malloc((numV + 1) * sizeof(int));		
		res = (complex<double>*)malloc(sizeof(complex<double>)*(int)(vartot));
		

		//Arreglo FOR que itera para los dos tipos de polarización
		for (int pol = 0; pol < 2; pol++)
		{
			cout << "Polarizacion: " << pol + 1 << endl;
			

			vector<vector<complex<double> > >KM(vartot);
			vector<vector<int> >ID(vartot);
			vector<vector<int> >ID2(vartot);
			vector<vector<int> >ID3(vartot);
			vector<vector<int> >DIR(vartot);
			vector<complex<double> >vecM;
			vector<doublecomplex>AR;
			vector<int>ja2;
			

			doublecomplex*aR;
			int*ja;

			/*PARTE 4: Enumeración de todas las aristas en el dominio*/
			q = (int)(eleY*nx*nz);
			q2 = q + (int)(eleZ*ny*nx);
			for (int k = 0; k < eleZ; k++)
			{
				for (int i = 0; i < eleX; i++)
				{
					for (int j = 0; j < eleY; j++)
					{
						p = (int)(k*(eleX*eleY) + i * (eleY)+j);
						lx[p] = (double)abs(x[i + 1] - x[i]);
						ly[p] = (double)abs(y[j + 1] - y[j]);
						lz[p] = (double)abs(z[k] - z[k + 1]);
						numb[p][0] = q2 + (int)(i + j * eleX + (k + 1)*eleX*ny) + 1;
						numb[p][1] = q2 + (int)(i + (j + 1) * eleX + (k + 1)*eleX*ny) + 1;
						numb[p][2] = q2 + (int)(i + j * eleX + (k)*eleX*ny) + 1;
						numb[p][3] = q2 + (int)(i + (j + 1) * eleX + (k)*eleX*ny) + 1;
						numb[p][4] = (int)(j + i * (eleY)+(k + 1)*(eleY*nx)) + 1;
						numb[p][5] = (int)(j + i * (eleY)+(k)*(eleY*nx)) + 1;
						numb[p][6] = (int)(j + (i + 1) * (eleY)+(k + 1)*(eleY*nx)) + 1;
						numb[p][7] = (int)(j + (i + 1) * (eleY)+(k)*(eleY*nx)) + 1;
						numb[p][8] = q + (int)(i*(eleZ*ny) + j * eleZ + k) + 1;
						numb[p][9] = q + (int)((i + 1)*(eleZ*ny) + j * eleZ + k) + 1;
						numb[p][10] = q + (int)((i)*(eleZ*ny) + (j + 1) * eleZ + k) + 1;
						numb[p][11] = q + (int)((i + 1)*(eleZ*ny) + (j + 1) * eleZ + k) + 1;

					}
				}
			}
			/*PARTE 5: Definición de los valores de resistividad en cada celda*/
			for (int k = 0; k < eleZ; k++)
			{
				for (int i = 0; i < eleX; i++)
				{
					for (int j = 0; j < eleY; j++)
					{

						p = (int)(j + i * eleY + (eleX*eleY)*k);
						if (k < num1z)
						{
							Resis[p] = ResA;
						}
						else
						{
							xmed = (x[i] + x[i + 1]) / 2.0;
							ymed = (y[j] + y[j + 1]) / 2.0;
							zmed = (z[k] + z[k + 1]) / 2.0;
							rdis = abs(sqrt(pow(xmed, 2.0) + pow(ymed, 2.0) + pow(zmed + 1.05, 2.0)));
							if ((CASO==1)&&(xmed >= -0.75) && (xmed <= 0.75)&&(ymed >= -0.75) && (ymed <= 0.75)&&(zmed<=-0.3)&&(zmed>=-1.8))
							{
								Resis[p] = 0.5;
							}
							else
							{
								Resis[p] = ResH;
							}
							if ((CASO==3)&&(rdis <= 0.75))
							{
								Resis[p] = 0.5;
							}
							else
							{
								Resis[p] = ResH;
							}
							if (CASO == 2)
							{	
								Resis[p] = ResH;
							}
							
							

						}
					}
				}
			}

			//Arriba dir Y
			for (int i = 0; i < nx; i++)
			{
				for (int j = 0; j < eleY; j++)
				{
					p = (int)(i*eleY) + j;
					if (pol == 0)
					{
						dirichV[p] = (complex<double>) A *(complex<double>) exp((complex<double>)k0*(complex<double>)z[0] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[0] * valC);
						cond[p] = 1;
					}
					else
					{
						dirichV[p] = 0.0;
						cond[p] = 1;
					}
				}
			}
			//Arriba dir X

			for (int j = 0; j < ny; j++)
			{
				for (int i = 0; i < eleX; i++)
				{
					p = (int)(eleY*nx*nz + eleZ * ny*nx) + (int)(i + eleX * j);
					if (pol == 1)
					{
						dirichV[p] = (complex<double>) A *(complex<double>) exp((complex<double>)k0*(complex<double>)z[0] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[0] * valC);
					}
					else
					{
						dirichV[p] = 0.0;
					}
					cond[p] = 1;

				}
			}
			//Izquierda dir Y
			for (int k = 1; k < (nz - 1); k++)
			{
				for (int j = 0; j < eleY; j++)
				{
					p = j + (int)(eleY*nx*k);
					if (k <= num1z)
					{
						if (pol == 0)
						{
							dirichV[p] = (complex<double>)A *(complex<double>) exp((complex<double>)k0*(complex<double>)z[k] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0* (complex<double>)z[k] * valC);
							//dirichI[q] = p + 1;
							cond[p] = 1;
						}
						else
						{
							dirichV[p] = 0.0;
							cond[p] = 1;
						}

					}
					else
					{
						if (pol == 0)
						{
							dirichV[p] = (complex<double>)C * (complex<double>)exp((complex<double>)k1*(complex<double>)z[k] * valC);
							//dirichI[q] = p + 1;
							cond[p] = 1;
						}
						else
						{
							dirichV[p] = 0.0;
							cond[p] = 1;
						}
					}
				}
			}
			//cout << "C: " << A << B << C << endl;
			//Derecha dir Y
			for (int k = 1; k < (nz - 1); k++)
			{
				for (int j = 0; j < eleY; j++)
				{
					p = j + (int)(eleY*nx*k + eleY * (nx - 1));

					if (k <= num1z)
					{
						if (pol == 0)
						{
							dirichV[p] = (complex<double>)A * (complex<double>)exp((complex<double>)k0*(complex<double>)z[k] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[k] * valC);
							//dirichI[q] = p + 1;
							cond[p] = 1;
						}
						else
						{
							dirichV[p] = 0.0;
							cond[p] = 1;
						}
					}
					else
					{
						//k0 = sqrt(im*mu0 / ResH * wr);
						if (pol == 0)
						{
							dirichV[p] = (complex<double>)C * (complex<double>)exp((complex<double>)k1*(complex<double>)z[k] * valC);
							cond[p] = 1;
						}
						else
						{
							dirichV[p] = 0.0;
							cond[p] = 1;
						}
					}

				}
			}


			//Atrás y adelante dir X
			for (int j = 0; j < 2; j++)
			{
				for (int k = 1; k < (nz - 1); k++)
				{
					for (int i = 0; i < eleX; i++)
					{
						p = (int)(eleY*nx*nz + eleZ * ny*nx) + (int)(i + eleX * ny* k + (ny - 1)*eleX*j);
						if (k <= num1z)
						{
							if (pol == 1)
							{
								dirichV[p] = (complex<double>)A * (complex<double>)exp((complex<double>)k0*(complex<double>)z[k] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[k] * valC);
								cond[p] = 1;
							}
							else
							{
								dirichV[p] = 0.0;
								cond[p] = 1;
							}
						}
						else
						{
							if (pol == 1)
							{
								dirichV[p] = (complex<double>)C * (complex<double>)exp((complex<double>)k1*(complex<double>)z[k] * valC);
								cond[p] = 1;
							}
							else
							{
								dirichV[p] = 0.0;
								cond[p] = 1;
							}
						}
					}
				}
			}	

			
			for (int i = 0; i < vartot; i++)
			{
				IA[i] = 0;
				JA[i] = 0;
			}
			IA[vartot] = 0;
			/*Parte 6: Ensamblaje no ordenado de las matrices elementales */
			enume = num1z * eleY*eleX;
			for (int p = 0; p < elemtot; p++)
			{
				for (int j = 0; j < 12; j++)
				{
					l[j] = (int)(numb[p][j] - 1);
				}
				for (int j = 0; j < 12; j++)
				{
					for (int i = 0; i < 12; i++)
					{
						if ((i < 4) && (j < 4))
						{
							valR = lx[p] * lz[p] / (6.0*ly[p])*K1[j][i] + lx[p] * ly[p] / (6.0*lz[p])*K2[j][i];
							valI = lx[p] * ly[p] * lz[p] / 36.0*F[j][i] * wr*mu0 / (Resis[p]);
						}
						else
						{
							if ((i >= 4) && (i < 8) && (j >= 4) && (j < 8))
							{
								valR = lx[p] * ly[p] / (6.0*lz[p])*K1[j - 4][i - 4] + ly[p] * lz[p] / (6.0*lx[p])*K2[j - 4][i - 4];
								valI = lx[p] * ly[p] * lz[p] / 36.0*F[j - 4][i - 4] * wr*mu0 / (Resis[p]);
							
							}
							else
							{
								if ((i >= 8) && (j >= 8))
								{
									valR = ly[p] * lz[p] / (6.0*lx[p])*K1[j - 8][i - 8] + lx[p] * lz[p] / (6.0*ly[p])*K2[j - 8][i - 8];
									valI = lx[p] * ly[p] * lz[p] / 36.0*F[j - 8][i - 8] * wr*mu0 / (Resis[p]);
									
								}
								else
								{
									valI = 0.0;
									if (((i >= 8) && (j >= 4)) || ((j >= 8) && (i >= 4)))
									{
										if (j < 8)
										{
											valP = K3[j - 4][i - 8];
										}
										else
										{
											valP = K3[i - 4][j - 8];
										}
										valR = -lx[p] / 6.0*valP;

									}
									if (((i >= 8) && (j < 4)) || ((i < 4) && (j >= 8)))
									{
										//
										if (i < 4)
										{
											valP = K3[j - 8][i];
										}
										else
										{
											valP = K3[i - 8][j];//j i-8
										}
										valR = -ly[p] / 6.0*valP;

									}
									if (((i < 8) && (j < 4)) || ((i < 4) && (j < 8)))
									{
										if (i >= 4)
										{
											valP = K3[j][i - 4];//j i-8
										}
										else
										{
											valP = K3[i][j - 4];
										}
										valR = -lz[p] / 6.0*valP;

									}
								}
							}

						}
						K[j][i].real((double)valR*(double)pow(10.0, 3.0));
						K[j][i].imag((double)valI*(double)pow(10.0, 9.0));
					}
				}
				for (int j = 0; j < 12; j++)
				{
					if (cond[l[j]] == 0)
					{
						for (int i = 0; i < 12; i++)
						{
							sizeF = ID[l[j]].size();
							valF = distance(ID[l[j]].begin(), find(ID[l[j]].begin(), ID[l[j]].end(), l[i] + 1));
							if ((sizeF == 0) || (valF == sizeF))
							{
								//cout << p << endl;
								KM[l[j]].push_back(K[j][i]);
								ID[l[j]].push_back(l[i] + 1);
								ID3[l[j]].push_back(ID[l[j]].size() - 1);
								if (cond[l[i]] != 0)
								{
									DIR[l[j]].push_back(l[i] + 1);
									ID2[l[j]].push_back(ID[l[j]].size());
								}
							}
							else
							{
								//cout << p<< endl;
								KM[l[j]][valF] = KM[l[j]][valF] + K[j][i];
								//ID[l[j]].push_back(l[i] + 1);
							}
						}
					}
				}
			}
			cout << "Done Assembling matrix" << endl;			
			ac = 0;
			id2 = 0;
			id = 0;
			max = 0;
			/*Parte 7: Construcción del vector B */
			for (int p = 0; p < vartot; p++)
			{
				DV[p] = 0.0;
				if (cond[p] == 0)
				{
					Sd = 0.0;
					if (DIR[p].size() > max)
					{
						max = DIR[p].size();
					}
					for (int i = 0; i < DIR[p].size(); i++)
					{
						id = DIR[p][i] - 1;
						id2 = ID2[p][i] - 1;
						//valF = distance(ID[p].begin(), find(ID[p].begin(), ID[p].end(), id + 1));
						Sd = Sd + (complex<double>) KM[p][id2] * (complex<double>)dirichV[id];// *0.0;// dirichV[id];
						KM[p][id2].real(0.0);
						KM[p][id2].imag(0.0);
					}
					DV[p] = DV[p] - Sd;
					varc[p] = ac;
				}
				else
				{
					ac = ac + 1;
					varc[p] = ac;
				}
			}

			cout << "Done correcting matrix" << endl;			
			numS = 0;
			//numS = 0;
			q = 0;
			p2 = -1;
			/*Parte 8: Reordenamiento de la matriz K */

			for (int p = 0; p < vartot; p++)
			{
				if (cond[p] == 0)
				{
					p2 = p2 + 1;
					vecID = ID[p];
					vecID3 = ID3[p];
					std::vector<int> vec(vecID3.size());
					IncGenerator g(0);
					std::generate(vec.begin(), vec.end(), g);
					//iota(vecID3.begin(), vecID3.end(), 0);
					sort(vec.begin(), vec.end(),myfunction);
					//	[&vecID](size_t i1, size_t i2) {return vecID[i1] < vecID[i2]; });
					for (int i = 0; i < ID[p].size(); i++)
					{
						id = vec[i];
						id2 = ID[p][id] - 1;
						//cout <<"vector"<< id << endl;
						if (id2 >= p)
						{
							if (abs(KM[p][id]) > pow(10.0, -10.0))
							{

								refV.re = KM[p][id].real();
								refV.i = KM[p][id].imag();
								AR.push_back(refV);
								ja2.push_back(id2 - varc[id2]);//ID[p][id]-1
								numS = numS + 1;
								q = q + 1;
							}
						}
						if ((id2) == p)
						{
							ia[p2] = (numS - 1);
						}
					}
					bv[p2].re = DV[p].real();
					bv[p2].i = DV[p].imag();
				}
			}

			cout << "Done reordering matrix" << endl;
			/*Parte 9: Solución del sistema de ecuaciones AX=B con PARDISO MKL*/
			cout << "Matrix dim: " << p2 + 1 << endl;
			ia[numV] = (numS);
			
			numT = q;
			
			aR = (doublecomplex*)malloc(numT * sizeof(doublecomplex));
			ja = (int*)malloc(numT * sizeof(int));
			for (int i = 0; i < numT; i++)
			{

				aR[i].re = AR[i].re;
				aR[i].i = AR[i].i;
				ja[i] = ja2[i];
			}
			int      nnz = ia[numV];
			int      mtype = 6;
			int      nrhs = 1;
			void    *pt[64];
			int      iparm[64];
			int      maxfct, mnum, phase, error, msglvl, solver;
			double   dparm[64];
			int      num_procs;
			doublecomplex   ddum;
			int	      idum;
			char *pValue;
			size_t len;
			char*var;
			int      i;	
			for (i = 0; i < 64; i++) {
				pt[i] = 0;
			}
			for (i = 0; i < 64; i++) {
				iparm[i] = 0;
			}  		
			mkl_set_num_threads(18);
			iparm[0]=1;
			iparm[1]=2;
			iparm[2] = 18;//(int)pValue;// num_procs;
			maxfct = 1;
			mnum = 1;
			msglvl = 1;
			error = 0;
			for (i = 0; i < numV + 1; i++) {
				ia[i] += 1;
			}
			for (i = 0; i < nnz; i++) {
				ja[i] += 1;
			}
			error = 0;
			solver = 0;
			//pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
			//pardiso_chkmatrix_z(&mtype, &numV, aR, ia, ja, &error);
			//pardiso_chkvec_z(&numV, &nrhs, bv, &error);
			//pardiso_printstats_z(&mtype, &numV, aR, ia, ja, &nrhs, bv, &error);
			phase = 11;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase,
				&numV, aR, ia, ja, &idum, &nrhs,
				iparm, &msglvl, &ddum, &ddum, &error);
			if (error != 0) {
				cout<<("\nERROR during symbolic factorization: %d", error)<<endl;
			exit(1);
			}
			phase = 22;
			cout<<("\nReordering completed ... ")<<endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase,
				&numV, aR, ia, ja, &idum, &nrhs,
				iparm, &msglvl, &ddum, &ddum, &error);
			if (error != 0) {
			cout<<("\nERROR during numerical factorization: %d", error)<<endl;
			exit(2);
			}
			cout << ("\nFactorization completed ...\n ") << endl;
			phase = 33;
			iparm[7] = 2;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase,
				&numV, aR, ia, ja, &idum, &nrhs,
				iparm, &msglvl, bv, xR, &error);
			if (error != 0) {
				cout<<("\nERROR during solution: %d", error)<<endl;
				exit(3);
			}
			cout << ("\nSolve completed ... ") << endl;
			//complex<double>refT(pow(10.0, 6.0), 0);


			p2 = 0;
			//ofstream myfilea;
			
			//myfilea.open("ResultadoModeloCirc5BB.txt");
			
			for (int i = 0; i < vartot; i++)
			{
				if (cond[i] == 0) {
					res[i].real(xR[p2].re);
					res[i].imag(xR[p2].i);
					p2 = p2 + 1;
				}
				else
				{
					res[i].real(dirichV[i].real());
					res[i].imag(dirichV[i].imag());
				}
				//myfilea << setprecision(10)<<res[i].real();
				//myfilea << " ";
				//myfilea << setprecision(10)<<res[i].imag();
				//myfilea << "\n";
			}

			//myfilea.close();
			/*ofstream myfile2a;
			myfile2a.open("EnumModeloCirc5BB.txt");
			myfile2a << eleX;
			myfile2a << "\n";
			myfile2a << eleY;
			myfile2a << "\n"; 
			myfile2a << eleZ;
			myfile2a << "\n";
			myfile2a << num1z;
			myfile2a << "\n";
			for (int i = 0; i < nx; i++)
			{
				myfile2a << x[i];
				myfile2a << "\n";
			}
			for (int i = 0; i < ny; i++)
			{
				myfile2a << y[i];
				myfile2a << "\n";
			}
			for (int i = 0; i < nz; i++)
			{
				myfile2a << z[i];
				myfile2a << "\n";
			}
			myfile2a.close();*/

			/*Parte 10: Calculo de resistividad aparente y fase en la superficie con diferencias finitas  */
			/*Parte 11: Generación de archivos de salida */
			int c = 0;
			ifstream myfile("ResistividadesCirc5VB.txt");
			ofstream myfile2, myfile3;
			if (pol == 1)
			{
				if (myfile)
				{
					myfile2.open("ResistividadesCirc5VB.txt", ios_base::app);
					myfile3.open("FasesCirc5VB.txt", ios_base::app);
				}
				else
				{
					myfile2.open("ResistividadesCirc5VB.txt");
					myfile3.open("FasesCirc5VB.txt");
				}
				//myfile2 << Freq;
				//myfile2 << " ";
				//myfile3 << Freq;
				//myfile3 << " ";
			}
			for (int i = 1; i < eleX; i++)
			{
				for (int j = 1; j < eleY; j++)
				{
					skY = (int)(j + i * eleY + num1z * (eleY*nx) - 1);
					skZ = (int)(eleY*nx*nz + num1z + i * (eleZ*ny) + j * (eleZ));
					skX = (int)(eleY*nx*nz + eleZ * ny*nx + num1z * eleX*ny + i + j * eleX - 1);
					a[0] = (double)(x[i] - x[i - 1])*1000.0;
					a[1] = (double)(y[j] - y[j - 1])*1000.0;
					b[0] = (double)(x[i + 1] - x[i])*1000.0;
					b[1] = (double)(y[j + 1] - y[j])*1000.0;
					for (int l = 0; l < 2; l++)
					{
						A1[l] = b[l] / (a[l] + b[l]);
						B1[l] = a[l] / (a[l] + b[l]);
						A2[l] = -1.0 / (a[l] * (a[l] / b[l] + 1.0));
						C2[l] = a[l] / (b[l] * b[l] * (a[l] / b[l] + 1.0));
						B2[l] = -A2[l] - C2[l];
					}
					if (pol == 0)
					{
						Ex1[c] = (A1[0] * res[skX] + B1[0] * res[skX + 1]);
						Ey1[c] = (A1[1] * res[skY] + B1[1] * res[skY + 1]);
						Ezy = (A2[1] * (res[skZ - eleZ - 1] + res[skZ - eleZ]) / 2.0 + B2[1] * (res[skZ] + res[skZ - 1]) / 2.0 + C2[1] * (res[skZ + eleZ - 1] + res[skZ + eleZ]) / 2.0);
						Ezx = (A2[0] * res[skZ - (int)(eleZ*ny)] + B2[0] * res[skZ] + C2[0] * res[skZ + (int)(eleZ*ny)]);
						Eyz = ((A1[1] * res[skY - (int)(eleY*nx)] + B1[1] * res[skY - (int)(eleY*nx) + 1]) - (A1[1] * res[skY + (int)(eleY*nx)] + B1[1] * res[skY + (int)(eleY*nx) + 1])) / ((double)(z[num1z - 1] - z[num1z + 1])*1000.0);
						Exz = ((A1[0] * res[skX] + B1[0] * res[skX + 1]) - (A1[0] * res[skX + (int)(eleX*ny)] + B1[0] * res[skX + (int)(eleX*ny + 1)])) / ((double)(z[num1z] - z[num1z + 1])*1000.0);
						//im = sqrt(-1.0);
						Hx1[c] = im * (complex<double>)(Eyz - Ezy) /(complex<double>) (wr*mu0);
						Hy1[c] = im * (complex<double>)(Exz - Ezx) / (complex<double>)(wr*mu0);
						
					}
					else
					{
						//im = sqrt(-1.0);
						Ex2[c] = (A1[0] * res[skX] + B1[0] * res[skX + 1]);
						Ey2[c] = (A1[1] * res[skY] + B1[1] * res[skY + 1]);
						Ezy = (A2[1] * (res[skZ - eleZ - 1] + res[skZ - eleZ]) / 2.0 + B2[1] * (res[skZ] + res[skZ - 1]) / 2.0 + C2[1] * (res[skZ + eleZ - 1] + res[skZ + eleZ]) / 2.0);
						Ezx = (A2[0] * (res[skZ - (int)(eleZ*ny) - 1] + res[skZ - (int)(eleZ*ny)]) / 2.0 + B2[0] * (res[skZ] + res[skZ - 1]) / 2.0 + C2[0] * (res[skZ + (int)(eleZ*ny) - 1] + res[skZ + (int)(eleZ*ny)]) / 2.0);

						Eyz = ((A1[1] * res[skY - (int)(eleY*nx)] + B1[1] * res[skY - (int)(eleY*nx) + 1]) - (A1[1] * res[skY + (int)(eleY*nx)] + B1[1] * res[skY + (int)(eleY*nx) + 1])) / ((double)(z[num1z - 1] - z[num1z + 1])*1000.0);
						Exz = ((A1[0] * res[skX - (int)(eleX*ny)] + B1[0] * res[skX - (int)(eleX*ny) + 1]) - (A1[0] * res[skX + (int)(eleX*ny)] + B1[0] * res[skX + (int)(eleX*ny + 1)])) / ((double)(z[num1z - 1] - z[num1z + 1])*1000.0);
						Hx2[c] = im * (Eyz - Ezy) / (wr*mu0);
						Hy2[c] = im * (Exz - Ezx) / (wr*mu0);
						
						detZ = (Hx2[c] * Hy1[c] - Hx1[c] * Hy2[c]);
						
						Zxx = (Ex2[c] * Hy1[c] - Ex1[c] * Hy2[c]) / detZ;
						Zxy = (-Ex2[c] * Hx1[c] + Ex1[c] * Hx2[c]) / detZ;
						Zyx = (Ey2[c] * Hy1[c] - Ey1[c] * Hy2[c]) / detZ;
						Zyy = (-Ey2[c] * Hx1[c] + Ey1[c] * Hx2[c]) / detZ;
						RhoXY = (double)pow(abs(Zxy), 2.0) / (wr*mu0);
						RhoYX = (double)pow(abs(Zyx), 2.0) / (wr*mu0);
						RhoYY = (double)pow(abs(Zyy), 2.0) / (wr*mu0);
						RhoXX = (double)pow(abs(Zxx), 2.0) / (wr*mu0);
						FaseXY = atan(imag(Zxy) / real(Zxy))*180.0 / M_PI;
						FaseYX = atan(imag(Zyx) / real(Zyx))*180.0 / M_PI;
						myfile2 << (RhoXY);
						myfile2 << " ";
						myfile2 << (RhoYX);
						myfile2 << " ";
						myfile2 << (RhoXX);
						myfile2 << " ";
						myfile2 << (RhoYY);
						myfile2 << "\n";
						myfile3 << FaseXY;
						myfile3 << " ";
						myfile3 << FaseYX;
						myfile3 << "\n";
						if ((i == (int)(eleX / 2)) && (j == (int)(eleY / 2)))
						{
							cout << "ResXY en centro " << RhoXY << endl;
							cout << "ResYX en centro " << RhoYX << endl;
							cout << "FaseXY en centro " << FaseXY << endl;
							cout << "FaseYX en centro " << FaseYX << endl;
						}
					}
					c = c + 1;
				}
			}
			if (pol == 1)
			{
				myfile2.close();
				myfile3.close();
			}
			/* -------------------------------------------------------------------- */
			/* ..  Convert matrix back to 0-based C-notation.                       */
			/* -------------------------------------------------------------------- */
			for (i = 0; i < numV + 1; i++) {
				ia[i] -= 1;
			}
			for (i = 0; i < nnz; i++) {
				ja[i] -= 1;
			}
			phase = -1;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase,
				&numV, &ddum, ia, ja, &idum, &nrhs,
				iparm, &msglvl, &ddum, &ddum, &error);

			KM.erase(KM.begin(), KM.end());
			ID.erase(ID.begin(), ID.end());
			ID2.erase(ID2.begin(), ID2.end());
			ID3.erase(ID3.begin(), ID3.end());
			DIR.erase(DIR.begin(), DIR.end());
			vecM.erase(vecM.begin(), vecM.end());
			AR.erase(AR.begin(), AR.end());
			ja2.erase(ja2.begin(), ja2.end());
			vecID.erase(vecID.begin(), vecID.end());
			vecID3.erase(vecID3.begin(), vecID3.end());
			free(aR);
			free(ja);
			
		}
		zP1.erase(zP1.begin(), zP1.end());
		zP2.erase(zP2.begin(), zP2.end());
		zP3.erase(zP3.begin(), zP3.end());
		free(dirichI);
		free(cond);
		free(dirichV);
		free(numb);
		free(lx);
		free(ly);
		free(lz);
		free(Resis);
		free(res);
		free(ia);
		free(varc);
		free(IA);
		free(JA);
		free(DV);
		free(xR);
		free(bv);
		free(z);
	
	}


	return 0;
}
