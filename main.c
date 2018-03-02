
#include <stdio.h>
#include <stdlib.h>

int InputFunction(char *inputfile, double ***q, int *nq, double **v, int dimension);
int nx1dInterpolation(double ** v, int * nq, double dq, int dimension, int n_spline);

int main(int argc, char **argv){

    int i, j, k;
    int dimension = 2;
    int n_points  = 0;
    int n_spline  = 1;

    if(argc > 2){
        dimension = atoi(argv[2]);
    }
    if(argc > 3){
        n_spline = atoi(argv[3]);
    }

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

// start interpolation process
    n_points = nx1dInterpolation(&v, nq, dq, dimension, n_spline);


// output
    printf("N");
    for(i = 0; i < dimension; ++i){ printf("\t% d", nq[i]); } printf("\n");

    switch(dimension){
        case 1:
        // 1D Output
            for(i = 0; i < n_points; ++i){
                printf("\t%d\t% lf\n", i+1, v[i]);
            }
            break;

        case 2:
        // 2D Output
            for(i = 0; i < nq[0]; ++i){
                for(j = 0; j < nq[1]; ++j){
                    printf("\t% 3d", i+1);
                    printf("\t% 3d", j+1);
                    printf("\t% .12lf", v[i*nq[1] + j]);
                    printf("\n");
                }
                printf("\n");
            }
            break;

        case 3:
        // 3D Output
            for(i = 0; i < nq[0]; ++i){
                for(j = 0; j < nq[1]; ++j){
                    for(k = 0; k < nq[2]; ++k){
                        printf("\t% 3d", i+1);
                        printf("\t% 3d", j+1);
                        printf("\t% 3d", k+1);

                        printf("\t% .12lf", v[k + nq[2]*j + nq[2]*nq[1]*i]);
                        printf("\n");
                    }
                    printf("\n");
                }
                printf("\n");
            }
            break;

        default:
            for(i = 0; i < n_points; ++i){
                printf("\t% .12lf\n", v[i]);
            }
            break;
    }

    return 0;
}
