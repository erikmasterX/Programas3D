
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
#include<numeric>
#include <stdlib.h>
#include <omp.h>
#include <iomanip>
#include <mkl_pardiso.h>
#include <mkl_types.h>
#include <mkl.h>
#define M_PI  3.14159265358979323846

using namespace std;
struct IncGenerator {
	int current_;
	IncGenerator(int start) : current_(start) {}
	int operator() () { return current_++; }
};
typedef struct {
	double re;
	double i;
}
doublecomplex;
vector<int>vecID, vecID3;
bool myfunction(int i, int j) { return (vecID[i] < vecID[j]); };
float Lx = 30.0, Ly = 30.0, Lz = 30.0; //Dimensiones del dominio en Km
int eleX, eleY, eleZ, num1z, num2z,num3z,num4z,num5z,q;
float *x, *y, *z, dx, dy, dz, z1, z2, x1, x2, dp;

int num1x, num2x, num3x, num4x, num5x;
int num1y, num2y, num3y, num4y, num5y;
int nx, ny, nz;
void CuerpoCircular();
void CuerpoHomogeneo();
void CuerpoRectangular();
void CuerpoCircular()
{
	/* Definición del numero de elementos y nodos en la malla */
	num1x = 10;
	num2x = 16;
	num1y = 10;
	num2y = 16;
	eleX = (int)(2 * num1x + num2x);
	eleY = (int)(2 * num1x + num2x);
	num1z = 15;
	num2z = 10;
	num3z = 16;
	num4z = 20;
	eleZ = (int)(num1z + num2z + num3z + num4z);
	cout << "eleX: " << eleX << "eleY: " << eleY << "eleZ: " << eleZ << endl;	
	nx = eleX + 1;
	ny = eleY + 1;
	nz = eleZ + 1;
	
	z = (float*)malloc(nz * sizeof(float));
	y = (float*)malloc(ny * sizeof(float));
	x = (float*)malloc(nx * sizeof(float));
	/* Definición de coordenadas en la dirección X */
	int p1;
	x1 = log10(1.0);
	x2 = log10(Lx / 2);
	dp = (x2 - x1) / (float)(num1x);
	p1 = 0;
	for (int i = 0; i < (num1x + 1); i++)
	{
		x[p1] = -pow(10.0, x2 - (double)i * dp);
		p1++;
	}
	dx = 2.0 / (float)(num2x);
	for (int i = 1; i < (num2x + 1); i++)
	{
		x[p1] = -1.0 + (float)dx*i;
		p1++;
	}

	for (int i = 1; i < (num1x + 1); i++)
	{
		x[p1] = pow(10.0, x1 + (double)i*dp);
		p1++;
	}
	/* Definición de coordenadas en la dirección Y */
	for (int i = 0; i < ny; i++)
	{
		y[i] = x[i];
	}
	/* Definición de coordenadas en la dirección Z */
	p1 = 0;
	z1 = log10(0.001);
	z2 = log10(10.0);
	dz = (z2 - z1) / (float)(num1z - 1);
	for (int i = 0; i < (num1z); i++)
	{
		z[p1] = (float)pow(10.0, z2 - (double)dz*(i));
		p1++;
	}
	z[p1] = 0;
	p1++;
	z1 = log10(0.001);
	z2 = log10(0.05);
	dz = (z2 - z1) / (float)(num2z - 1);
	for (int i = 0; i < num2z; i++) {
		z[p1] = -(float)pow(10.0, z1 + (double)dz*i);
		p1++;
	}
	dz = 2.0 / (float)num3z;
	for (int i = 1; i < (num3z + 1); i++) {
		z[p1] = -(0.05 + (float)dz*i);
		p1++;
	}
	z1 = log10(2.05);
	z2 = log10(15.0);
	dz = (z2 - z1) / (float)num4z;
	for (int i = 1; i < (num4z + 1); i++) {
		z[p1] = -pow(10.0, z1 + (double)dz*i);
		p1++;
	}

}
void CuerpoRectangular() {

	/*Definición del número de elementos y nodos en el dominio*/
	num1x = 10;
	num2x = 8;
	num3x = 10;
	num1y = 10;
	num2y = 8;
	num3y = 10;
	eleX = (int)(2 * num1x + 2 * num2x + num3x);
	eleY = (int)(2 * num1x + 2 * num2x + num3x);
	num1z = 25;//38
	num2z = 10;
	num3z = 8;
	num4z = 10; 
	num5z = 25;//num5z=14;
	eleZ = (int)(num1z + num2z + 2 * num3z + num4z + num5z);
	cout << "eleX: " << eleX << "eleY: " << eleY << "eleZ: " << eleZ << endl;
	nx = eleX + 1;
	ny = eleY + 1;
	nz = eleZ + 1;
	
	z = (float*)malloc(nz * sizeof(float));
	y = (float*)malloc(ny * sizeof(float));
	x = (float*)malloc(nx * sizeof(float));

	/*Definición de coordenadas en X*/
	int p1;
	x1 = log10(0.95);
	x2 = log10(Lx / 2);
	dp = (x2 - x1) / (float)(num1x);
	p1 = 0;
	for (int i = 0; i < (num1x + 1); i++)
	{
		x[p1] = -pow(10.0, x2 - (double)i * dp);
		p1++;
	}
	dx = 0.4 / (float)(num2x);
	for (int i = 1; i < (num2x + 1); i++)
	{
		x[p1] = -0.95 + (float)dx*i;
		p1++;
	}
	dx = (0.55 + 0.55) / (float)num3x;
	for (int i = 1; i < (num3x + 1); i++)
	{
		x[p1] = -0.55 + (float)dx*i;
		p1++;
	}
	dx = 0.4 / (float)num2x;
	for (int i = 1; i < (num2x + 1); i++)
	{
		x[p1] = 0.55 + (float)dx*i;
		p1++;
	}
	for (int i = 1; i < (num1x + 1); i++)
	{
		x[p1] = pow(10.0, x1 + (double)i*dp);
		p1++;
	}
	/*Definición de coordenadas en Y*/
	for (int i = 0; i < ny; i++)
	{
		y[i] = x[i];
	}

	/*Definición de coordenadas en Z*/
	p1 = 0;
	z1 = log10(0.001);
	z2 = log10(10.0);
	dz = (z2 - z1) / (float)(num1z - 1);
	for (int i = 0; i < (num1z); i++)
	{
		z[p1] = (float)pow(10.0, z2 - (double)dz*(i));
		p1++;
	}
	z[p1] = 0;
	p1++;
	z1 = log10(0.001);
	z2 = log10(0.1);
	dz = (z2 - z1) / (float)(num2z - 1);
	for (int i = 0; i < num2z; i++) {
		z[p1] = -(float)pow(10.0, z1 + (double)dz*i);
		p1++;
	}
	dz = 0.4 / (float)num3z;
	for (int i = 1; i < (num3z + 1); i++) {
		z[p1] = -(0.1 + (float)dz*i);
		p1++;
	}
	dz = 1.1 / num4z;
	for (int i = 1; i < (num4z + 1); i++) {
		z[p1] = -(0.5 + (float)dz*i);
		p1++;
	}
	dz = 0.4 / num3z;
	for (int i = 1; i < (num3z + 1); i++) {
		z[p1] = -(1.6 + (float)dz*i);
		p1++;
	}
	z1 = log10(2);
	z2 = log10(Lz / 2);
	dz = (z2 - z1) / (float)num5z;
	for (int i = 1; i < (num5z + 1); i++) {
		z[p1] = -pow(10.0, z1 + (double)dz*i);
		p1++;
	}

}
void CuerpoHomogeneo()
{
	/*Definición  del número de elementos y nodos en el dominio*/
	eleX = 35;
	eleY = 35;
	num1z = 35;
	num2z = 35;
	eleZ = num1z + num2z;

	nx = eleX + 1;
	ny = eleY + 1;
	nz = eleZ + 1;
	
	dx = Lx / (float)(nx - 1);
	dy = Ly / (float)(ny - 1);

	z = (float*)malloc(nz * sizeof(float));
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
	/*Definición de coordenadas en Z*/
	q = 0;
	z1 = log10(0.001);
	z2 = log10(Lz / 2.0);
	dz = (z2 - z1) / (float)(num1z - 1);
	for (int i = 0; i < (num1z); i++)
	{
		z[q] = (float)pow(10.0, z2 - (double)dz*(i));
		q = q + 1;
	}
	z[q] = 0;
	q = q + 1;
	dz = (z2 - z1) / (float)(num2z - 1);
	for (int i = 0; i < num2z; i++) {
		z[q] = -(float)pow(10.0, z1 + (double)dz*i);
		q = q + 1;
	}

}
int CASO;
int main()
{	
	int aristasFronNum, elemtot, vartot;
	int *dirichI, p = 0,*cond;
	complex <float> *dirichV;
	complex<float> H0(1.0, 0.0), A, B, C;
	double *Resis, ResA;
	double xmed, ymed, zmed, rdis;

	/* PARTE 1: Creación de malla rectilínea para cualquiera de 3 casos:
	   - CASO=1: Prisma rectangular embebido en medio homogéneo
	   - CASO=2: Medio homogéneo
	   - CASO=3: Cuerpo "esférico" embebido en medio homogéneo */ 
	CuerpoRectangular(); //Caso de prisma
	CASO = 1;

	aristasFronNum = (int)((eleY * nz + eleZ * ny) * 2 +((nx - 2)*eleZ+nz*eleX) * 2 + (eleY*(nx - 2) + eleX * (ny- 2)) * 2); //Número de aristas en la frontera
	elemtot = (int)(eleX*eleY*eleZ); //Número de prismas en todo el dominio
	vartot = (int)(((eleY*nz) + (eleZ)*ny)*nx+eleX*ny*nz);//Número de aristas en todo el dominio
	//dirichI = (int*)malloc(aristasFronNum * sizeof(int)); 
	cond = (int*)malloc(vartot * sizeof(int)); //Arreglo indicador de fronteras
	double Freq = 10.0, wr, ResH = 100.0; //Frecuencia y resistividad de medio.
	//cout << elemtot << endl;
	Resis = (double*)malloc(elemtot * sizeof(double)); //Arreglo con valores de resistividad de cada celda	
	ResA = (double)3.3*pow(10.0, 8.0); //Resisitividad del Aire

	/*PARTE 2: Definición de la resistividad de cada celda*/
	for (int k = 0; k < eleZ; k++)
	{
		for (int i = 0; i < eleX; i++)
		{
			for (int j = 0; j < eleY; j++)
			{				
				p = (int)(j + i * eleY + (eleX*eleY)*k);
				if (k < num1z)
				{
					//Definición de la resistividad en el aire
					Resis[p] = ResA;
				}
				else
				{
					//Definición de la resistividad en el subsuelo
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
						Resis[p] = 0.5;
					}
			        if ((CASO==3)&&(rdis<=0.75))
					{
						Resis[p] = 0.5;
					}					
					else
					{
						Resis[p] = ResH;
					}	
					if ((CASO == 2))
					{
						Resis[p] = ResH;
					}
				}
			}
		}
	}

	/*PARTE 3: Definición de los valores en las fronteras*/
	/*Solución analítica aplicada en las fronteras*/
	complex<float>k0, k1, val1, im;	
	complex<double>valC(1000.0, 0.0);	
	double mu0;
	double valY;
	im = sqrt(-1);
	mu0 = 4.0 * M_PI*pow(10.0, -7.0);
	wr = 2.0 * M_PI*Freq;
	complex<float>vecR1(0.0, mu0 / ResA * wr), vecR2(0.0, mu0 / ResH * wr);
	complex<float>val2(1000.0, 0.0);
	complex<float>varA(-1.0, 0.0);
	k0 = sqrt(vecR1);
	k1 = sqrt(vecR2);
	C = H0;	
	A.real((double)(H0.real()/2.0+H0.real()*pow(10.0,4.0)*sqrt(3.3/ResH)/2.0));	
	A.imag(0.0);
	B.real(C.real() - A.real());
	B.imag(C.imag()-A.imag());

	dirichV = (complex<float>*)malloc(vartot * sizeof(complex<float>));
	q = 0;
	
	for (int i = 0; i < vartot; i++)
	{
		cond[i] = 0;
	}
	//Arriba dir Y
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < eleY; j++)
		{
			p = (int)(i*eleY) + j;
			dirichV[p] =0.0;// (complex<double>) A *(complex<double>) exp((complex<double>)k0*(complex<double>)z[0] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[0] * valC);
			cond[p] = 1;
			q = q + 1;
		}
	}
	//Arriba dir X
	
	for (int j = 0; j < ny; j++)
	{
		for (int i = 0; i < eleX; i++)
		{
			p = (int)(eleY*nx*nz+eleZ*ny*nx)+(int)(i+eleX*j);
			dirichV[p] =(complex<double>) A *(complex<double>) exp((complex<double>)k0*(complex<double>)z[0] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[0] * valC);
			cond[p] = 1;

		}
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
			p = (int)(nz-1)*(eleY*nx) + (int)(i*eleY) + j;
			dirichV[p] = 0.0;
			cond[p] = 1;
	
			q = q + 1;
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
				dirichV[p] = 0.0;// (complex<double>)A *(complex<double>) exp((complex<double>)k0*(complex<double>)z[k] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0* (complex<double>)z[k] * valC);
				//dirichI[q] = p + 1;
				cond[p] = 1;
			}
			else
			{
				dirichV[p] =0.0; //(complex<double>)C * (complex<double>)exp((complex<double>)k1*(complex<double>)z[k] * valC);
				//dirichI[q] = p + 1;
				cond[p] = 1;
			}
		}
	}	
	//Derecha dir Y
	for (int k = 1; k < (nz - 1); k++)
	{
		for (int j = 0; j < eleY; j++)
		{
			p = j + (int)(eleY*nx*k + eleY * (nx - 1));
			
				if (k <= num1z)
				{

					dirichV[p] = 0.0; //(complex<double>)A * (complex<double>)exp((complex<double>)k0*(complex<double>)z[k] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[k] * valC);
					//dirichI[q] = p + 1;
					cond[p] = 1;					
				}
				else
				{
					//k0 = sqrt(im*mu0 / ResH * wr);
					dirichV[p] = 0.0;//(complex<double>)C * (complex<double>)exp((complex<double>)k1*(complex<double>)z[k] * valC);
					//dirichI[q] = p + 1;					
					cond[p] = 1;
				}
		
		}
	}
	
	complex<float>val3(0.0, 0.0);
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
	//Atrás y adelante dir X
	for (int j = 0; j < 2; j++)
	{
		for (int k = 1; k < (nz - 1); k++)
		{
			for (int i = 0; i < eleX; i++)
			{
				p = (int)(eleY*nx*nz + eleZ * ny*nx) + (int)(i + eleX *ny* k + (ny - 1)*eleX*j);
				if (k <= num1z)
				{
					dirichV[p] = //(complex<double>)A * (complex<double>)exp((complex<double>)k0*(complex<double>)z[k] * valC) + (complex<double>)B *(complex<double>) exp(-(complex<double>)k0 * (complex<double>)z[k] * valC);
					cond[p] = 1;
				}
				else
				{
					dirichV[p] =0.0; //(complex<double>)C * (complex<double>)exp((complex<double>)k1*(complex<double>)z[k] * valC);
					cond[p] = 1;
				}
				
				//cond[p] = 1;
			}
		}
	}
	//int chk;
	//chk = (int)(eleY*nx*nz+eleZ*ny*nx+10);

	
	cout << "VariablesNum:" << vartot << endl;
	cout << "FronterasNum:" << aristasFronNum << endl;
	
	int **numb;
	double*lx, *ly, *lz;
	lx = (double*)malloc(elemtot * sizeof(double));
	ly = (double*)malloc(elemtot * sizeof(double));
	lz = (double*)malloc(elemtot * sizeof(double));

	numb = (int**)malloc(elemtot * sizeof(int *));
	int q2;

	/*PARTE 4: Enumeración de las aristas en todo el dominio*/

	for (int i = 0; i < elemtot; i++)
	{
		numb[i] = (int*)malloc(12 * sizeof(int));
	}
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
				numb[p][0] = q2 + (int)(i + j * eleX + (k + 1)*eleX*ny)+1;
				numb[p][1] = q2 + (int)(i + (j + 1) * eleX + (k + 1)*eleX*ny)+1;
				numb[p][2] = q2 + (int)(i + j * eleX + (k)*eleX*ny)+1;
				numb[p][3] = q2 + (int)(i + (j + 1) * eleX + (k)*eleX*ny)+1;
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
	double F[4][4] = {  {4.0,2.0,2.0,1.0},
					   {2.0,4.0,1.0,2.0},
					   {2.0,1.0,4.0,2.0},
					   {1.0,2.0,2.0,4.0} };
	//cout << K3[1][0] << endl;
	complex<double>K[12][12];
	double valP,valR,valI;
	vector<vector<complex<double> > >KM(vartot);
	vector<vector<int> >ID(vartot);
	vector<vector<int> >ID2(vartot);
	vector<vector<int> >ID3(vartot);
	vector<vector<int> >DIR(vartot);
	
	int *indX, *indY,*IA,*JA;
	indX = (int*)malloc(vartot * sizeof(int));
	indY= (int*)malloc(vartot* sizeof(int));
	IA= (int*)malloc((vartot+1) * sizeof(int));
	JA = (int*)malloc(vartot * sizeof(int));
	
	
	for (int i = 0; i < vartot; i++)
	{
		IA[i] = 0;
		JA[i] = 0;
		//indX[i] = 0;
		//indY[i] = 0;
	}
	IA[vartot] = 0;
	
	int valF, sizeF;
	int enume;
	enume = num1z * eleY*eleX;
	/*Parte 5: Ensamblaje no ordenado de las matrices elementales */
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
					valI = lx[p] * ly[p] * lz[p] / 36.0*F[j][i] * wr*mu0 /(Resis[p]);
					/*if (p < enume)
					{
						valI = 0.0;
					}*/
					
				}
				else
				{
					if ((i >= 4) && (i < 8) && (j >= 4) && (j < 8))
					{
						valR = lx[p] * ly[p] / (6.0*lz[p])*K1[j-4][i-4] + ly[p] * lz[p] / (6.0*lx[p])*K2[j-4][i-4];
						valI = lx[p] * ly[p] * lz[p] / 36.0*F[j-4][i-4] * wr*mu0 / (Resis[p]);
						/*if (p < enume)
						{
							valI = 0.0;
						}*/
					}
					else
					{
						if ((i >= 8) && (j >= 8))
						{
							valR = ly[p] * lz[p] / (6.0*lx[p])*K1[j - 8][i - 8] + lx[p] * lz[p] / (6.0*ly[p])*K2[j - 8][i - 8];
							valI = lx[p] * ly[p] * lz[p] / 36.0*F[j - 8][i - 8] * wr*mu0 / (Resis[p]);
							/*if (p < enume)
							{
								valI = 0.0;
							}*/
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
							if (((i >= 8) && (j <4)) || ((i < 4) && (j >= 8)))
							{
								//
								if (i<4)
								{
									valP = K3[j-8][i];
								}
								else
								{
									valP = K3[i-8][j];//j i-8
								}
								valR = -ly[p] / 6.0*valP;

							}
							if (((i < 8) && (j < 4)) || ((i < 4) && (j < 8)))
							{
								if (i >=4)
								{
									valP = K3[j][i-4];//j i-8
								}
								else
								{
									valP = K3[i][j-4];
								}
								valR = -lz[p] / 6.0*valP;

							}
						}
					}

				}	
				
				K[j][i].real((double)valR*(double)pow(10.0,3.0));
				K[j][i].imag((double)valI*(double)pow(10.0,9.0));
				
			}
			
		}
		
		for (int j = 0; j < 12; j++)
		{
			if (cond[l[j]] == 0)
			{
				for (int i = 0; i < 12; i++)
				{
					sizeF = ID[l[j]].size();
					valF = distance(ID[l[j]].begin(),find(ID[l[j]].begin(),ID[l[j]].end(),l[i]+1));
					if ((sizeF==0)||(valF==sizeF))
					{
						//cout << p << endl;
						KM[l[j]].push_back(K[j][i]);
						ID[l[j]].push_back(l[i] + 1);
						ID3[l[j]].push_back(ID[l[j]].size()-1);
						if (cond[l[i]] != 0)
						{
							DIR[l[j]].push_back(l[i] + 1);
							ID2[l[j]].push_back(ID[l[j]].size());
						}
					}
					else
					{
						//cout << p<< endl;
						KM[l[j]][valF] = KM[l[j]][valF]+K[j][i];
						//ID[l[j]].push_back(l[i] + 1);
					}
				}
			}
		}
	}
	/*Parte 6: Construcción del vector B */

	cout << "Ensamblado terminado" << endl;
	complex<double>Sd;
	complex<double>*DV;
	DV = (complex<double>*)malloc(vartot*sizeof(complex<double>));
	int id = 0,id2=0,max,ac=0;
	max = 0;
	int*varc;
	varc = (int*)malloc(vartot*sizeof(int));

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
				
				
				Sd = Sd + (complex<double>) KM[p][id2]*(complex<double>)dirichV[id];
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
	
	doublecomplex*xR,*bv;
	int*ia;
	int    numV;
	numV = (int)(vartot - aristasFronNum);
	cout << numV << endl;
	xR = (doublecomplex*)malloc(numV * sizeof(doublecomplex));
	bv = (doublecomplex*)malloc(numV * sizeof(doublecomplex));
	ia = (int*)malloc((numV+1) * sizeof(int));
	
	vector<complex<double> >vecM;
	vector<doublecomplex>AR;
	vector<int>ja2;
	
	doublecomplex refV;
	int numS = 0;
	//numS = 0;
	q = 0;
	int p2 = -1;

	/*Parte 7: Reordenación de la matriz K */
	
	for (int p = 0; p < vartot; p++)
	{
		if (cond[p] == 0)
		{
			p2 = p2 + 1;
			vecID = ID[p];
			vecID3 = ID3[p];
			std::vector<int> vec(vecID3.size()); // vector with 100 ints.
			IncGenerator g(0);
			std::generate(vec.begin(), vec.end(), g);
			//std::iota(vecID3.begin(),vecID3.end(), 0);
			sort(vec.begin(), vec.end(), myfunction);//[&vecID](size_t i1, size_t i2) {return vecID[i1] < vecID[i2]; });
			
			for (int i = 0; i < ID[p].size(); i++)
			{
				id = vec[i];
				id2 = ID[p][id] - 1;
			
				if (id2 >= p)
				{
					if (abs(KM[p][id])>pow(10.0,-10.0))
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
					ia[p2]=(numS-1);					
				}
							
			}			
			bv[p2].re = DV[p].real();
			bv[p2].i = DV[p].imag();
		}		
	}
	cout << "Reordenamiento terminado " << endl;

	/*Parte 8: Solución del sistema de ecuaciones AX=B con PARDISO MKL*/
	
	ia[numV] =(numS);
	int numT;
	numT = q;
	doublecomplex*aR;
	int*ja;
	aR = (doublecomplex*)malloc(numT * sizeof(doublecomplex));
	ja = (int*)malloc(numT * sizeof(int));
	for (int i = 0; i < numT; i++)
	{		
		aR[i].re =  AR[i].re;
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
	int      i;
	doublecomplex   ddum;        /* Double dummy */
	int	      idum;              /* Integer dummy. */
	char *pValue;
	size_t len;
	char*var;
	maxfct = 1;        
	mnum = 1;   
	msglvl = 1;       
	error = 0;  
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}      
	iparm[0]=1;
	iparm[1]=2;

	mkl_set_num_threads(18); //Aquí se definie el número de hilos a usar.
	
	iparm[2] = 18;
	iparm[59]=0;

	for (i = 0; i < numV + 1; i++) {
		ia[i] += 1;
	}
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}
	
	error = 0;
	solver = 0; 	
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
	phase = 11;
	error=0;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&numV, aR, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		cout<<("\nERROR during symbolic factorization: %d", error)<<endl;
		exit(1);
	}
	cout<<("\nReordering completed ... ")<<endl;
	cout<<("\nNumber of nonzeros in factors  = %d", iparm[17])<<endl;
	cout<<("\nNumber of factorization MFLOPS = %d", iparm[18])<<endl;
	phase = 22;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&numV, aR, ia, ja, &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		cout<<("\nERROR during numerical factorization: %d", error)<<endl;
		exit(2);
	}
	cout<<("\nFactorization completed ...\n ")<<endl;
	phase = 33;
	iparm[7] =2;//1;      

	pardiso(pt, &maxfct, &mnum, &mtype, &phase,
		&numV, aR, ia, ja, &idum, &nrhs,
		iparm, &msglvl, bv, xR, &error);

	if (error != 0) {
		cout<<("\nERROR during solution: %d", error)<<endl;
		exit(3);
	}

	cout<<("\nSolve completed ... ")<<endl;

	complex<double>refT(pow(10.0, 6.0), 0);
	
	/*Parte 9: Escritura de archivos de salida */

	ofstream myfile;
	doublecomplex resul;

	/*Componentes reales e imaginarias en todos el espacio y en las 3 direcciones (X,Y,Z)*/
	myfile.open("ResultadoCamposE.txt");
	p2 = 0;
	for (int i = 0; i < vartot; i++)
	{
		if (cond[i] == 0) {
			resul.re = xR[p2].re;
			resul.i = xR[p2].i;
			p2 = p2 + 1;
		}
		else
		{
			resul.re = dirichV[i].real();
			resul.i = dirichV[i].imag();
		}
		myfile << setprecision(10)<<resul.re;
		myfile << " ";
		myfile << setprecision(10)<<resul.i;
		myfile << "\n";
	}
	
	myfile.close();
	
	/*Información respecto a los modelos usados*/
	ofstream myfile2;
	myfile2.open("InformacionModelo.txt");
	myfile2 << eleX;
	myfile2 << "\n";
	myfile2 << eleY;
	myfile2 << "\n"; 
	myfile2 << eleZ;
	myfile2 << "\n";
	myfile2 << num1z;
	myfile2 << "\n";
	for (int i = 0; i < nx; i++)
	{
		myfile2 << x[i];
		myfile2 << "\n";
	}
	for (int i = 0; i < ny; i++)
	{
		myfile2 << y[i];
		myfile2 << "\n";
	}
	for (int i = 0; i < nz; i++)
	{
		myfile2 << z[i];
		myfile2 << "\n";
	}
	myfile2.close();


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
	cout << "Terminado" << endl;
	
	return 0;
}


