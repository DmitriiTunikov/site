#pragma once
#define _USE_MATH_DEFINES 
#include <math.h> 

enum 
{
  even,
  cheb,
  my_grid
};

typedef struct grid_t
{
  int type;
  int pointsCount;
  double *x;
  double *y;
}grid_t;

typedef struct coef_t
{
  double a;
  double b;
  double c;
  double d;
}coef_t;

double efunc(double x);
double GetDeltaByCubicSplain(grid_t grid, double(*f)(double), int condition);
void GetGrid(grid_t *grid, int grid_type, double(*f)(double), double A, double B);
double GetDeltaByLagrangh(grid_t grid, double(*f)(double));