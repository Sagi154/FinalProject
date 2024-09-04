# include <stdio.h>
# include <stdlib.h>
# include <math.h>

struct cluster
{
 int vect_count;
 int vect_length;
 double* centroid;
 double* sum;
} ;

struct cord
{
 double value;
 struct cord *next;
};

struct vector
{
 struct vector *next;
 struct cord *cords;
};

int create_vector_arr(struct vector* head_vec);
double **symnmf();
/*
Perform full the symNMF as described in 1 and output H.
 */
double **sym();
/*
Calculate and output the similarity matrix as described in 1.1.
 */
double **ddg();
/*
Calculate and output the Diagonal Degree Matrix as described in 1.2
 */
double **norm();
/*
Calculate and output the normalized similarity matrix as described in 1.3
 */
