
#include <stdio.h>
#include <stdlib.h>

int InputFunction(char *inputfile, double ***q, int *nq, double **v, int dimension);
int nx1dInterpolation(double ** v, int * nq, double dq, int dimension, int n_spline);

int main(int argc, char **argv){

    int i, j, k, l;
    int jump;

/* Input */

    int dimension = 2;
    int n_points  = 0;
    int n_spline  = 1;

    double ** q = NULL;
    double  * v = NULL;
    char    * inputfile = argv[1];
    int     * nq = NULL;
    double    dq;

    if(argc > 2){
        dimension = atoi(argv[2]);
    }
    if(argc > 3){
        n_spline = atoi(argv[3]);
    }

// Memory allocation
    nq = calloc(dimension, sizeof(int));
    v  = malloc(sizeof(double));
    q  = malloc(sizeof(double*));
    for(i = 0; i < dimension; ++i){
        q[i] = malloc(sizeof(double));
    }

// Actual input
    n_points = InputFunction(inputfile, &q, nq, &v, dimension);

// check input
    for(i = 0, k = 1; i < dimension; ++i){
        k *= nq[i];
    }
    if(k != n_points){
        fprintf(stderr, "\n (-) Error in input file:"
                        "\n     Product of dimension lengths (%d) does not match n_points (%d)"
                        "\n     Aborting...\n\n"
                        , k, n_points
               );
        exit(1);
    }

// calculate dq, assume equi distant spacing
    dq = q[dimension-1][1] - q[dimension-1][0];


// start interpolation process
    n_points = nx1dInterpolation(&v, nq, dq, dimension, n_spline);

// calculate new dq
    dq = dq / (double)(n_spline + 1);

// reallocate memory for all q[*]s
    for(i = 0; i < dimension; ++i){
        q[i] = realloc(q[i], n_points * sizeof(double));
    }

// fill all q[*]s
    for(i = dimension-1, jump = 1; i >= 0; --i){

        for(j = 0; j < n_points/jump/nq[i]; ++j){
            for(k = 0; k < nq[i]; ++k){
                for(l = 0; l < jump; ++l){

                    q[i][l + k*jump + j*nq[i]*jump] = q[i][0] + (double)k * dq;

                }
            }
        }
        jump *= nq[i];
    }


// output
    printf("N");
    for(i = 0; i < dimension; ++i){ printf("\t% d", nq[i]); } printf("\n");

    for(i = 0, j = 0; i < n_points; ++i){
    // add blank lines
        for(j = (dimension-1), k = 1; j >= 0; --j){
            k *= nq[j];
            if(i%k == 0){ printf("\n"); }
        }
    // actual output
        for(j = 0; j < dimension; ++j){
            printf("\t% lf", q[j][i]);
        }
        printf("\t% lf", v[i]);
        printf("\n");
    }

    return 0;
}
