
#include <stdio.h>
#include <stdlib.h>

int InputFunction(char *inputfile, double ***q, int *nq, double **V, int dimension);
int interpolator(double ** v, int * nq, double dq, int dimension, int n_spline);

int main(int argc, char **argv){

    int i, j;
    int dimension = 1;
    int n_points  = 0;
    int n_spline  = 1;

/* Input */

    double ** q = NULL;
    double  * v = NULL;
    char    * inputfile = argv[1];
    int     * nq = NULL;
    double dq;

// Memory allocation
    nq = calloc(dimension, sizeof(int));
    v  = malloc(sizeof(double));
    q  = malloc(sizeof(double*));
    for(i = 0; i < dimension; ++i){
        q[i] = malloc(sizeof(double));
    }

// Actual input
    n_points = InputFunction(inputfile, &q, nq, &v, dimension);
    dq = q[dimension-1][1] - q[dimension-1][0];


/*
// output input
    for(i = 0; i < n_points; ++i){
        printf("\t%d\t% lf\n", i, v[i]);
    }
    printf("\n");
//*/


// start interpolation process
    n_points = interpolator(&v, nq, dq, dimension, n_spline);

    for(i = 0; i < n_points; ++i){
        printf("\t%d\t% lf\n", i+1, v[i]);
    }

    return 0;
}


#include <stdio.h>
#include <stdlib.h>
int interpolator(double ** v, int * nq, double dq, int dimension, int n_spline){

/*  One dimensional cubic spline interpolation:
{{{
    The 1D cubic spline interpolation is based on the calculation of a third
    order (cubic) polynomial between each subset of neighbouring grid data points.

        f(h(i))   =   a  +   b h(i)  +   c h(i)^2  + d h(i)^3
        f'(h(i))  =   b  + 2 c h(i)  + 3 d h(i)^2
        f''(h(i)) = 2 c  + 6 d h(i)

    Where h is the grid spacing (h_n = x_n+1 - x_n). Using the first and second
    derivatives and the natural boundary conditions f''(h_0) = 0 and f''(h_n) = 0
    allows for the derivation of the following recursion formulas:

        a(i+1) = a(i) +   b(i) h(i) +   c(i) h(i)^2 + d(i) h(i)^3
        b(i+1) = b(i) + 2 c(i) h(i) + 3 d(i) h(i)^2
        c(i+1) = c(i) + 3 d(i) h(i)

    Here a directly corresponds to the original data points, giving 3n equations for
    the calculation of 3n unknown parameters. Rearrangement leads to

        3 ( (a(i+1) - a(i)) / h(i) - (a(i) - a(i-1)) / h(i-1)
        =
        h(i-1) c(i-1) + 2 (h(i) + h(i-1)) c(i) + h(i) c(i+1) 

    On an equispaced grid (as required by the Numerov method) this simplifies to

        3/h^2 ( a(i-1) - 2 a(i) + a(i+1) ) = c(i-1) + 4 c(i) + c(i+1)

    In matrix notation (Ax = b) this can be represented as tridiagonal matrix equation

        | 4  1            | |  c_1  |          |             0              |
        | 1  4  1         | |  c_2  |          | a(0)   - 2 a(1)   + a(2)   |
        |    .  .  .      | |   .   |  = 3/h/h |             .              |
        |      .  .  .    | |   .   |          |             .              |
        |         1  4  1 | | c_n-2 |          | a(n-3) - 2 a(n-2) + a(n-1) |
        |            1  4 | | c_n-1 |          | a(n-2) - 2 a(n-1) + a(n)   |

    which can be solved using the Thomas (tridiagonal matrix) algorithm (vide infra).

    With this information the b and d coefficients can be calculated

        d(i) = 1/(3 h) * (c(i+1) - c(i))
        b(i) = 1/h     * (a(i+1) - a(i)) - h/3 * (2 c(i) + c(i+1))

//}}}*/

/*  Thomas algorithm for solving tridiagonal matrix problems:
 {{{
    The Thomas algorithm is a very efficient way to solve equations of the type
 
        a_i x_i-1 + b_i x_i + c_i x_i+1 = d_i
 
    with the two constraints a_1 = 0 and c_n = 0
 
    In matrix notation (Ax = b) this corresponds to
 
       | b_1    c_1                        | | x_1 |     | d_1 |
       | a_2    b_2    c_2                 | | x_2 |     | d_2 |
       |         .      .      .           | |  .  |  =  |  .  |
       |                .      .      .    | |  .  |     |  .  |
       |                      a_n    b_n   | | x_n |     | d_n |
 
    The algorithm starts with a forward sweep to eliminate the
    a coefficients by modifying the coefficients c (to c') and
    d (to d')
 
        c'_{i} =          c_{i}           / (b_{i} - a_{i} c'_{i-1})
        d'_{i} = (d_{i} - a_{i} d'_{i-1}) / (b_{i} - a_{i} c'_{i-1})
 
    Note, since a_1 = 0:
        c'_1 = c_1 / b_1
        d'_1 = d_1 / b_1
 
    The solution is then obtained by back substitution
 
        x_{i} = d'_{i} - c'_{i} x_{i+1}
 
    Since c_n = 0:
        x_{n} = d'_{n}
 
//}}}*/

/*
   In this particular case
       ∀ i ∊ (1,n] : a_i = 1
       ∀ i ∊ [1,n] : b_i = 4
       ∀ i ∊ [1,n) : c_i = 1

   allowing for further simplifications
     i = 1
       c'_{1} = 1/4
       d'_{1} = 0

     i ∊ [2,n)
       c'_{i} = 1 / (4 - c'_{i-1})
       d'_{i} = (d_{i} - d'_{i-1}) * c'_{i}

     i = n
       c'_{n} = 0
       d'_{n} = (d_{n} - d'_{n-1}) / (4 - c'_{n-1})
//*/

//----------------------------------------------------------------------
//  Declaration   Declaration   Declaration   Declaration   Declaration
//----------------------------------------------------------------------
// For better readability the parameters of the Thomas algorithm will
//  be prefixed with "matrix_" and the interpolation parameters will
//  be prefixed with "inter_"

    int i, j;
    int maxdim   = 0;
    int n_points = 0;
    int nn_points = 0;
    double newdq = 0;

// tridiagonal matrix algorithm coefficients
    double * matrix_c = NULL;
    double * matrix_d = NULL;
// interpolation coefficients
    double * inter_c  = NULL;
    double   inter_b;
    double   inter_d;
// interpolated data array
    double * yy = NULL;


//----------------------------------------------------------------------
//  Initialisation   Initialisation   Initialisation   Initialisation
//----------------------------------------------------------------------
// determine the size of the dimension with the most entries
//  and calculate n_points
    for(i = 0, maxdim = 0, n_points = 1, nn_points = 1; i < dimension; ++i){
        if(nq[i] > maxdim){
            maxdim = nq[i];
        }
        n_points *= nq[i];
        nn_points *= ((nq[i] - 1) * (n_spline + 1) + 1);
    }

// allocate memory for matrix_c, matrix_d and inter_c arrays
//  as well as the new array containing interpolated data yy
    matrix_c = malloc((maxdim-1) * sizeof(double));
    matrix_d = malloc((maxdim-1) * sizeof(double));
    inter_c  = malloc((maxdim-1) * sizeof(double));
    yy       = malloc(nn_points  * sizeof(double));

// define new delta q (newdq) which accounts for interpolated points
    newdq = dq / ((double)n_spline + 1.0);


//----------------------------------------------------------------------
//   Forward sweep   Forward sweep   Forward sweep   Forward sweep
//----------------------------------------------------------------------
// fill matrix_d array with tridiagonal matrix d elements
    matrix_d[0] = 0.0;
    for(i = 1; i < nq[0]-1; ++i){
        matrix_d[i] = 3/dq/dq * ((*v)[i-1] - 2*(*v)[i] + (*v)[i+1]);
    }

// calculate tridiagonal matrix elements c' and d' and save them to
//  matrix_c and matrix_d, respectively (overwrite matrix_d)
    matrix_c[0] = 1.0/4.0;
    for(i = 1; i < nq[0]-2; ++i){
        matrix_c[i] = 1.0 / (4.0 - matrix_c[i-1]);
        matrix_d[i] = (matrix_d[i] - matrix_d[i-1]) * matrix_c[i];
    }


//----------------------------------------------------------------------
// Back substitution to determine inter_c coefficients and interpolation
//----------------------------------------------------------------------
// solve the last entry of the matrix problem
    inter_c[i] = (matrix_d[i] - matrix_d[i-1]) / (4.0 - matrix_c[i-1]);

// calculate inter_b and inter_d values from inter_c
    inter_d = - inter_c[i] / (3.0 * dq);
    inter_b = ((*v)[i+1] - (*v)[i]) / dq - (2.0 / 3.0 * dq * inter_c[i]);

// interpolation procedure for the last n_spline+1 points
    for(j = n_spline+1; j >= 0; --j){
        yy[i * (n_spline+1) + j] = (*v)[i]
                                 + inter_b    * (newdq*j)
                                 + inter_c[i] * (newdq*j) * (newdq*j)
                                 + inter_d    * (newdq*j) * (newdq*j) * (newdq*j);
    }


    for(i = nq[0]-3; i >= 0; --i){
    // solve the rest of the matrix problem from back to front
        inter_c[i] = matrix_d[i] - matrix_c[i] * inter_c[i+1];

    // calculate inter_b and inter_d values from inter_c
        inter_d = (inter_c[i+1] - inter_c[i]) / (3.0 * dq);
        inter_b = ((*v)[i+1] - (*v)[i]) / dq - dq / 3.0 * (2.0 * inter_c[i] + inter_c[i+1]);


    // interpolation procedure for the remaining points (all but the last set)
        for(j = n_spline; j >= 0; --j){
            yy[i * (n_spline+1) + j] = (*v)[i]
                                     + inter_b    * (newdq*j)
                                     + inter_c[i] * (newdq*j) * (newdq*j)
                                     + inter_d    * (newdq*j) * (newdq*j) * (newdq*j);
        }
    }


// free original v array and point it to yy instead
    free(*v); *v = NULL;
    *v = yy;

// free memory
    free(matrix_c); matrix_c = NULL;
    free(matrix_d); matrix_d = NULL;
    free(inter_c);  inter_c  = NULL;

printf("#dimensions:");
for(i = 0; i < dimension; ++i){ printf("\t%d", nq[i]); } printf("\n");
printf("#maxdim    = \t%d\n", maxdim);
printf("#   dq     = \t% lf\n", dq);
printf("#newdq     = \t% lf\n", newdq);
printf("#n_spline  = \t%d\n", n_spline);
printf("#n_points  = \t%d\n", n_points);
printf("#nn_points = \t%d\n", nn_points);

    return nn_points;

}
