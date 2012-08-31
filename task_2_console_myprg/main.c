//
//  main.c
//  task2_console
//
//  Created by Alexander Demidov on 12.04.12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N1 32
#define N2 50
#define J 100
#define SHOW_K 99 //временной слой для вывода
#define SHOW_N1 5 //слой по иксу для вывода
#define SHOW_N2 5 //слой по игреку для вывода

#define X_MAX 2.0
#define X_MIN 0.0
#define Y_MAX M_PI
#define Y_MIN 0.0

#define T_MAX 0.003
#define T_MIN 0.0

//граница по иксу
#define CAPPA_X 1.0
#define NU_X 0.0
//граница по игреку
#define CAPPA_Y 1.0
#define NU_Y 0.0

#define SPEED 1.0

//начальные условия
double G(double x, double y)
{
    return cos(3.0*M_PI*x)*cos(2.0*y);
}

//нелинейная часть уравнения
double F(double t, double x, double y)
{
    return 0.0;
}

void solveMatrix(const int n, double A, double B, double C, double *f, double *x, const double cappa, const double nu, double *alfa, double *betta)
{
    alfa[0] = betta[0] = 0.0;
    alfa[1] = cappa;
    betta[1] = nu;
    int i;
    for (i = 1; i < (n-1); i++) {
        alfa[i+1] = B / (C - A*alfa[i]);
        betta[i+1] = (A*betta[i] + f[i]) / (C - A*alfa[i]);
    }
    x[n-1] = (nu + cappa*betta[n-1]) / (1.0 - cappa*alfa[n-1]);
    for (i = n-2; i >= 0; i--) {
        x[i] = alfa[i+1]*x[i+1] + betta[i+1];
    }
}

int main(int argc, const char * argv[])
{
	const double dX = (X_MAX - X_MIN) / N1;
	const double dY = (Y_MAX - Y_MIN) / N2;
	const double dT = (T_MAX - T_MIN) / J;
	double currentX;
	double currentY;
    double currentT;
	int i, j, k;
    //сеточная функция
	double ***GRID;
	GRID = (double ***)malloc(N1 * sizeof(double **));
    if (GRID == NULL) {
        printf("Memory error.\n");
        return 2;
    }
	for (i = 0; i < N1; i++) {
		GRID[i] = (double **)malloc(N2 * sizeof(double *));
        if (GRID[i] == NULL) {
            printf("Memory error.\n");
            return 2;
        }
		for (j = 0; j < N2; j++) {
			GRID[i][j] = (double *)malloc(J * sizeof(double));
            if (GRID[i][j] == NULL) {
                printf("Memory error.\n");
                return 2;
            }
		}
	}
    //задаём начальные условия
    //в начальный момент времени - функция G(x,y)
    currentX = X_MIN;
	for (i = 0; i < N1; i++) {
        currentY = Y_MIN;
		for (j = 0; j < N2; j++) {
			GRID[i][j][0] = G(currentX, currentY);
			currentY += dY;
		}
		currentX += dX;
	}
	
    //нужно сначала перейти от массива GRID[i][j][k] к temp[i][j] (переход от слоя k к k+1/2)
    //выделяем память для слоя k+1/2
	double **temp;
	temp = (double **)malloc(N1 * sizeof(double *));
    if (temp == NULL) {
        printf("Memory error.\n");
        return 2;
    }
	for (i = 0; i < N1; i++) {
		temp[i] = (double *)malloc(N2 * sizeof(double));
        if (temp[i] == NULL) {
            printf("Memory error.\n");
            return 2;
        }
	}
    
    //начинаем вычисления
	const double gamma1 = dT / (dX*dX);
	const double gamma2 = dT / (dY*dY);
    printf("gamma1 = %lf gamma2 = %lf\ndT = %lf dX*dX/4 = %lf dY*dY/4 = %lf dX = %lf dY = %lf\n", gamma1, gamma2, dT, dX*dX/4.0, dY*dY/4.0, dX, dY);
    
    double *fToTempLayer = (double *)malloc(sizeof(double) * N1); //столбец свободных членов при переходе на промежуточный слой
    double *xToTempLayer = (double *)malloc(sizeof(double) * N1); //неизвестный столбец при переходе на промежуточный слой
    //вспомогательные числа для прямой прогонки
    double *alfaToTempLayer = (double *)malloc(sizeof(double) * N1);
    double *bettaToTempLayer = (double *)malloc(sizeof(double) * N1);
    
    //всё аналогично для перехода с промежуточного на постоянный слой
    double *fFromTempLayer = (double *)malloc(sizeof(double) * N2);
    double *xFromTempLayer = (double *)malloc(sizeof(double) * N2);
    double *alfaFromTempLayer = (double *)malloc(sizeof(double) * N2);
    double *bettaFromTempLayer = (double *)malloc(sizeof(double) * N2);
    
    if (fToTempLayer == NULL || xToTempLayer == NULL || alfaToTempLayer == NULL || bettaToTempLayer == NULL || fFromTempLayer == NULL || xFromTempLayer == NULL || alfaFromTempLayer == NULL || bettaFromTempLayer == NULL) {
        printf("Memory error.\n");
        return 2;
    }
    //коэффициенты уравнения
    double AtoTemp, BtoTemp, CtoTemp;
    double AfromTemp, BfromTemp, CfromTemp;
    AtoTemp = BtoTemp = 0.5*gamma1*SPEED;
    CtoTemp = -(1.0 + gamma1*SPEED);
    AfromTemp = BfromTemp = 0.5*gamma2*SPEED;
    CfromTemp = -(1.0 + gamma2*SPEED);
    
    currentT = T_MIN;
    for (k = 1; k < J; k++) {
        //переходим на промежуточный слой
        currentY = Y_MIN + dY;
        for (j = 1; j < (N2-1); j++) {
            //заполняем столбец F
            currentX = X_MIN;
            for (i = 0; i < N1; i++) {
                fToTempLayer[i] = -(0.5*gamma2*SPEED*(GRID[i][j-1][k-1] + GRID[i][j+1][k-1]) + (1.0 - gamma2*SPEED)*GRID[i][j][k-1] + 0.5*dT*F(currentT + 0.5*dT, currentX, currentY));
                currentX += dX;
            }
            //решаем СЛАУ методом прогонки
            solveMatrix(N1, AtoTemp, BtoTemp, CtoTemp, fToTempLayer, xToTempLayer, CAPPA_X, NU_X, alfaToTempLayer, bettaToTempLayer);
            for (i = 0; i < N1; i++)
                temp[i][j] = xToTempLayer[i];
            currentY += dY;
        }
        for (i = 0; i < N1; i++) {
            temp[i][0] = CAPPA_Y*temp[i][1] + NU_Y;
            temp[i][N2-1] = CAPPA_Y*temp[i][N2-2] + NU_Y;
        }
        
        //переходим на слой k+1
        currentX = X_MIN + dX;
        for (i = 1; i < (N1-1); i++) {
            //заполняем столбец F
            currentY = Y_MIN;
            for (j = 0; j < N2; j++) {
                fFromTempLayer[j] = -(0.5*gamma1*SPEED*(temp[i-1][j] + temp[i+1][j]) + (1.0 - gamma1*SPEED)*temp[i][j] + 0.5*dT*F(currentT + 0.5*dT, currentX, currentY));
                currentY += dY;
            }
            //решаем СЛАУ методом прогонки
            solveMatrix(N2, AfromTemp, BfromTemp, CfromTemp, fFromTempLayer, xFromTempLayer, CAPPA_Y, NU_Y, alfaFromTempLayer, bettaFromTempLayer);
            for (j = 0; j < N2; j++) {
                GRID[i][j][k] = xFromTempLayer[j];
            }
            currentX += dX;
        }
        for (j = 0; j < N2; j++) {
            GRID[0][j][k] = CAPPA_X*GRID[1][j][k] + NU_X;
            GRID[N1-1][j][k] = CAPPA_X*GRID[N1-2][j][k] + NU_X;
        }
        currentT += dT;
    }
    
    //вывод в файл
    FILE *dataFile, *dataFile2, *dataFile3, *dataFile4;
    dataFile = fopen("points.txt", "w");
    dataFile2 = fopen("points2.txt", "w");
    dataFile3 = fopen("points3.txt", "w");
    dataFile4 = fopen("points4.txt", "w");
    if (dataFile == NULL || dataFile2 == NULL || dataFile3 == NULL || dataFile4 == NULL) {
        printf("Files openning error!\n");
        return 2;
    }
    
    currentY = Y_MIN;
    for(j = 0; j < N2; j++) {
        currentX = X_MIN;
		for(i = 0; i < N1; i++) {
			fprintf(dataFile2,"%.4lf\t%.3lf\t%.3lf\n", currentX, currentY, GRID[i][j][SHOW_K]);
            fprintf(dataFile,"%.4lf\t%.3lf\t%.3lf\n", currentX, currentY, GRID[i][j][0]);
            currentX += dX;
		}
        currentY += dY;
	}
    
    currentX = X_MIN;
    for (i = 0; i < N1; i++) {
        currentT = T_MIN;
        for (k = 0; k < J; k++) {
            fprintf(dataFile3,"%.4lf\t%.3lf\t%.3lf\n", currentX, currentT, GRID[i][SHOW_N2][k]);
            currentT += dT;
        }
        currentX += dX;
    }
    
    currentY = Y_MIN;
    for (j = 0; j < N2; j++) {
        currentT = T_MIN;
        for (k = 0; k < J; k++) {
            fprintf(dataFile4,"%.4lf\t%.3lf\t%.3lf\n", currentY, currentT, GRID[SHOW_N1][j][k]);
            currentT += dT;
        }
        currentY += dY;
    }
    
    fclose(dataFile);
    fclose(dataFile2);
    fclose(dataFile3);
    fclose(dataFile4);
    printf("Files are ready!\n");
    
    //освобождаем память
	for (i = 0; i < N1; i++) {
		for (j = 0; j < N2; j++) {
			free(GRID[i][j]);
			GRID[i][j] = NULL;
		}
		free(temp[i]);
		free(GRID[i]);
		GRID[i] = NULL;
		temp[i] = NULL;
	}
    free(fToTempLayer);
    free(xToTempLayer);
    free(alfaToTempLayer);
    free(bettaToTempLayer);
    free(fFromTempLayer);
    free(xFromTempLayer);
    free(alfaFromTempLayer);
    free(bettaFromTempLayer);
	free(GRID);
	free(temp);
    fToTempLayer = NULL;
    xToTempLayer = NULL;
    alfaToTempLayer = bettaToTempLayer = NULL;
    fFromTempLayer = NULL;
    xFromTempLayer = NULL;
    alfaFromTempLayer = bettaFromTempLayer = NULL;
	GRID = NULL;
	temp = NULL;
    return 0;
}
