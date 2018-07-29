#include <stdio.h>
#include "interpol.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

double fr(double x)
{
  return (1. / (25 + x * x));
}

int main()
{
  grid_t *grid;
  double delta = 0;
  char *filename = malloc(sizeof(char) * 20);
  int i = 0, points_count = 0, grids_count = 3, j;
  char *grid_type = malloc(sizeof(char) * 30);
  char *cur_f = malloc(sizeof(char) * 10);
  double(*f)(double) = NULL;
  //LabA
  for (j = 0; j < 2; j++)
  {
    //accept func
    if (j == 0)
    {
      strcpy(cur_f, "e^x");
      f = &efunc;
    }
    else 
    {
      strcpy(cur_f, "fr");
      f = &fr;
    }//get delta for diffrent grids
    for (i = 0; i < grids_count; i++)
    {
      switch (i)
      {
      case 0:
        strcpy(grid_type, "even");
        break;
      case 1:
        strcpy(grid_type, "chebish");
        break;
      case 2:
        strcpy(grid_type, "my_grid");
        break;
      }

      for (points_count = 10; points_count < 21; points_count += 4)
      {
        grid = malloc(sizeof(grid_t));
        grid->x = malloc(sizeof(double) * points_count);
        grid->y = malloc(sizeof(double) * points_count);
        grid->pointsCount = points_count;
        grid->type = i;
        if (f == &efunc)
          GetGrid(grid, i, f, -.92, 0.92);
        else
          GetGrid(grid, i, f, 0, M_PI * 0.3);
        delta = GetDeltaByLagrangh(*grid, f);
        printf("Function: %s   Grid type: %s   Point count: %d  Delta:%g\n", cur_f, grid_type, points_count, delta);
        free(grid->x);
        free(grid->y);
        free(grid);
      }
    }
  }

  printf("\nDelta expirements:\nFunction: tg   Grid type: even   Point count: 5\n");
  //зависимость от близости к pi/2 для tg
  for (double b_end = 9/10.0 * M_PI_2; b_end <= 9.9/10.0 * M_PI_2; b_end += 0.02 * M_PI_2)
  {
    grid = malloc(sizeof(grid_t));
    grid->x = malloc(sizeof(double) * 20);
    grid->y = malloc(sizeof(double) * 20);
    grid->pointsCount = 5;
    
    GetGrid(grid, 0, &tan, 0, b_end); 
    delta = GetDeltaByLagrangh(*grid, &tan);
    
    printf("delta:%g || Polinom delta:%g\n", M_PI_2 - b_end, delta);
    
    free(grid->x);
    free(grid->y);
    free(grid);
  }

  for (int cond = 0; cond < 2; cond++)
  {
  for (j = 0; j < 2; j++)
  {
    //accept func
    if (j == 0)
    {
      strcpy(cur_f, "Fr");
      f = &fr;
    }
    else
    {
      strcpy(cur_f, "tg");
      f = &tan;
    }//get delta for diffrent grids
    
    if (cond)
      printf("Not a knot condition:\n");
    for (points_count = 10; points_count < 16; points_count += 3)
    {
      grid = malloc(sizeof(grid_t));
      grid->x = malloc(sizeof(double) * points_count);
      grid->y = malloc(sizeof(double) * points_count);
      grid->pointsCount = points_count;
      grid->type = i;
      if (f == &fr)
        GetGrid(grid, 0, f, -1, 1);
      else
        GetGrid(grid, 0, f, 0, M_PI * 0.3);
      delta = GetDeltaByCubicSplain(*grid, f, cond);
      printf("Function: %s   Point count: %d  Delta:%g\n", cur_f, points_count, delta);
      free(grid->x);
      free(grid->y);
      free(grid);
    }
  }
  }
  
  //зависимость от близости к pi/2 для tg
  f = &tan;
  for (i = 5; i < 16; i+=5)
  {
    printf("points count: %d\n", i);
  for (double b_end = 9 / 10.0 * M_PI_2; b_end <= 9.9 / 10.0 * M_PI_2; b_end += 0.02 * M_PI_2)
  {
    grid = malloc(sizeof(grid_t));
    grid->x = malloc(sizeof(double) * 20);
    grid->y = malloc(sizeof(double) * 20);
    grid->pointsCount = i;
    GetGrid(grid, 0, &tan, 0, b_end);
    delta = GetDeltaByCubicSplain(*grid, f, 0);
  
    printf("delta:%g || Polinom delta:%g\n", M_PI_2 - b_end, delta);
  
    free(grid->x);
    free(grid->y);
    free(grid);
  }
  }
  return 0;
}