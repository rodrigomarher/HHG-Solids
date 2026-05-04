#include "vec3_util.h"

double dot(double *v1, double *v2){
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void cross(double *v1, double *v2, double *result){
    result[0] =  v1[1]*v2[2] - v1[2]*v2[1]; 
    result[1] = -v1[0]*v2[2] + v1[2]*v2[0]; 
    result[2] =  v1[0]*v2[1] - v1[1]*v2[0];
}

void solve_3by3(double **mat, double *v, double *result){
    double m00 = mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2];
    double m01 = mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2];
    double m02 = mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2];
    
    double detA =  mat[0][0]*m00
                 - mat[1][0]*m01
                 + mat[2][0]*m02;
    double detX =   v[0]*m00
                  - mat[1][0]*(v[1]*mat[2][2] - mat[2][1]*v[2])
                  + mat[2][1]*(v[1]*mat[1][2] - mat[1][1]*v[2]);
    double detY =   mat[0][0]*(v[1]*mat[2][2] - mat[2][1]*v[2])
                  - v[0]*m01
                  + mat[2][0]*(mat[0][1]*v[2] - v[1]*mat[0][2]);

    double detZ =  mat[0][0]*(mat[1][1]*v[2] - v[1]*mat[1][2])
                 - mat[1][0]*(mat[0][1]*v[2] - v[1]*mat[0][2])
                 + v[0]*m02;

    double invDet = 1.0 / detA;
    result[0] = detX*invDet;
    result[1] = detY*invDet;
    result[2] = detZ*invDet;
}
