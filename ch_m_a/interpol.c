#include "interpol.h"
#include <stdlib.h>
#include <stdio.h>

double efunc(double x)
{
  return exp(x) / (1 - x * x);
}

void GetGrid(grid_t *grid, int grid_type, double (*f)(double), double A, double B)
{
  int i;
  if (grid_type == even)
  {
    grid->x[0] = A;
    grid->y[0] = f(grid->x[0]);
    for (i = 1; i < grid->pointsCount; i++)
    {
      grid->x[i] = (A + B) * i / (grid->pointsCount - 1.0);
      grid->y[i] = f(grid->x[i]);
    }
    grid->x[grid->pointsCount - 1] = B;
    grid->y[grid->pointsCount - 1] = f(grid->x[grid->pointsCount - 1]);
  }
  else if (grid_type == cheb)
  {
    for (i = 1; i <= grid->pointsCount; i++)
    {
      grid->x[i-1] = 0.5 * ((A + B) + (B - A) * cos((2 * i - 1) / 2.0 / grid->pointsCount * M_PI));
      grid->y[i-1] = f(grid->x[i-1]);
    }
  }
  else if (grid_type == my_grid)
  {
    grid->x[0] = 0;
    grid->y[0] = f(grid->x[0]);
    for (i = 1; i < grid->pointsCount; i++)
    {
      if (f == &cos)
        grid->x[i] = A + cos(i-1) * M_PI_2;
      else
        grid->x[i] = A + cos(i -1) * M_PI_2 * 3/10.0;
      grid->y[i] = f(grid->x[i]);
    }
  }
}

double GetLagranghValue(double cur_x, grid_t grid)
{
  int i, j;
  double lagr_sum = 0, mult = 1;

  for (i = 0; i < grid.pointsCount; i++)
  {
    mult = 1;
    for (j = 0; j < grid.pointsCount; j++)
    {
      if (j != i)
        mult *= (cur_x - grid.x[j]) / (grid.x[i] - grid.x[j]);
    }
    lagr_sum += mult * grid.y[i];
  }
  return lagr_sum;
}

double fr(double x)
{
  return (1. / (25 + x * x));
}


double GetDeltaByLagrangh(grid_t grid, double (*f)(double))
{
  double delta = 0, new_delta = 0;
  double cur_x, cur_y;
  int i, k;
  FILE *fl = NULL, *fl_y = NULL;
  //Save data to file
  if (f == &fr && (grid.pointsCount == 10) && grid.type == even)
  {
    fl = fopen("frx.txt", "w");
    fl_y = fopen("fry.txt", "w");
    k = 0;
    for (cur_x = -.92; cur_x <= 0.92; cur_x += 1.84 / 100)
    {
      k++;
      cur_y = GetLagranghValue(cur_x, grid);
      fprintf(fl, "%lf ", cur_x);
      fprintf(fl_y, "%lf ", cur_y);
    }
  }
  else if (f == &tan && (grid.pointsCount == 3) && grid.type == even)
  {
    fl = fopen("tgX.txt", "w");
    fl_y = fopen("tgY3.txt", "w");
    for (cur_x = 0; cur_x <= M_PI_2; cur_x += M_PI_2 / 100.0)
    {
      cur_y = GetLagranghValue(cur_x, grid);
      fprintf(fl, "%lf ", cur_x);
      fprintf(fl_y, "%lf ", cur_y);
    }
  }

  for (i = 0; i < grid.pointsCount - 1; i++)
  {
    cur_x = (grid.x[i] + grid.x[i + 1]) / 2.0;
    cur_y = GetLagranghValue(cur_x, grid);
    new_delta = cur_y - f(cur_x);
    if (fabs(new_delta) > fabs(delta))
      delta = fabs(new_delta);
  }
  if (fl != NULL)
  {
   fclose(fl);
   fclose(fl_y);
  }
  return delta;
}

double GetSplainValue(coef_t *coefs, double x, grid_t grid)
{
  int i = 0, j = grid.pointsCount - 1;
  
  while (i + 1 < j)
  {
    int k = i + (j - i) / 2;
    if (x <= grid.x[k])
      j = k;
    else
      i = k;
  }
  
  double dx = (x - grid.x[j]);
  return coefs[j].a + (coefs[j].b + (coefs[j].c / 2. + coefs[j].d * dx / 6.) * dx) * dx; // Вычисляем значение сплайна в заданной точке.
}

double GetDeltaByCubicSplain(grid_t grid, double(*f)(double), int condition)
{
  int i, k;
  double delta1, delta2;
  FILE *fl, *fl_y;
  coef_t *coefs = malloc(sizeof(coef_t) * grid.pointsCount);
  double cur_x, cur_y, delta = 0, new_delta = 0;
  double *alpha = malloc(sizeof(double) * (grid.pointsCount - 1));
  double *beta = malloc(sizeof(double) * (grid.pointsCount - 1));
  double A, B, C, F, h_i, h_i1, z;

  //get a[i] coefs
  for (i = 0; i < grid.pointsCount; i++)
    coefs[i].a = f(grid.x[i]);
  if (!condition)
    coefs[0].c = 0;

  // Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
  // Вычисление прогоночных коэффициентов - прямой ход метода прогонки
  alpha[0] = beta[0] = 0.;
  int n = grid.pointsCount;

  for (i = 1; i < n - 1; ++i)
  {
    h_i = grid.x[i] - grid.x[i - 1], h_i1 = grid.x[i + 1] - grid.x[i];
    if (condition && (i == 1 || i == n - 2))//not a knot condition
    {
      if (i == 1)
      {
        A = 0;
        C = 3 * h_i + 2 * h_i1 + pow(h_i, 2)/ h_i1;
        B = h_i1;
      }
      else
      {
        A = h_i - pow(h_i1, 2)/h_i;
        C = 3 * h_i1 + 2 * h_i + pow(h_i1, 2) / h_i;
        B = 0;
      }
    }
    else
    {
      A = h_i;
      C = 2. * (h_i + h_i1);
      B = h_i1;
    }
    F = 6. * ((grid.y[i + 1] - grid.y[i]) / h_i1 - (grid.y[i] - grid.y[i - 1]) / h_i);
    z = (A * alpha[i - 1] + C);
    alpha[i] = -B / z;
    beta[i] = (F - A * beta[i - 1]) / z;
  }
  coefs[n - 1].c = (F - A * beta[n - 2]) / (C + A * alpha[n - 2]);

  // Нахождение решения - обратный ход метода прогонки
  for (int i = n - 2; i > 0; --i)
    coefs[i].c = alpha[i] * coefs[i + 1].c + beta[i];

  if (condition)
  {
    delta1 = grid.x[1] - grid.x[0];
    delta2 = grid.x[2] - grid.x[1];
    coefs[0].c = coefs[1].c * (3 * (delta1)+2 * delta2 + pow(delta1, 2) / delta2) + coefs[2].c * (delta2 - pow(delta1, 2) / delta2);
  }
  // По известным коэффициентам c[i] находим значения b[i] и d[i]
  for (int i = n - 1; i > 0; --i)
  {
    double h_i = grid.x[i] - grid.x[i - 1];
    
    coefs[i].d = (coefs[i].c - coefs[i - 1].c) / h_i;
    coefs[i].b = h_i * (2 * coefs[i].c + coefs[i - 1].c) / 6. + (grid.y[i] - grid.y[i - 1]) / h_i;
  }

  if (f == &fr && (grid.pointsCount == 10) && condition)
  {
    fl = fopen("frx10.txt", "w");
    fl_y = fopen("fry10.txt", "w");
    k = 0;
    for (cur_x = 0; cur_x <= M_PI_2; cur_x += M_PI_2 / 100.0)
    {
      k++;
      cur_y = GetSplainValue(coefs, cur_x, grid);
      fprintf(fl, "%lf ", cur_x);
      fprintf(fl_y, "%lf ", cur_y);
    }
  }
  else if (f == &tan && (grid.pointsCount == 5) && condition)
  {
    fl = fopen("tx5.txt", "w");
    fl_y = fopen("ty5.txt", "w");
    for (cur_x = 0; cur_x <= M_PI * 3 / 10.; cur_x += M_PI_2 / 100.0)
    {
      cur_y = GetSplainValue(coefs, cur_x, grid);
      fprintf(fl, "%lf ", cur_x);
      fprintf(fl_y, "%lf ", cur_y);
    }
  }

  //find delta
  for (i = 0; i < grid.pointsCount - 1; i++)
  {
    cur_x = (grid.x[i] + grid.x[i + 1]) / 2;
    cur_y = GetSplainValue(coefs, cur_x, grid);
    new_delta = cur_y - f(cur_x);
    if (fabs(new_delta) > fabs(delta))
      delta = fabs(new_delta);
  }

  free(coefs);
  free(alpha);
  free(beta);
  return delta;
}