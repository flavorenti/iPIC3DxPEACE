/* iPIC3D was originally developed by Stefano Markidis and Giovanni Lapenta. 
 * This release was contributed by Alec Johnson and Ivy Bo Peng.
 * Publications that use results from iPIC3D need to properly cite  
 * 'S. Markidis, G. Lapenta, and Rizwan-uddin. "Multi-scale simulations of 
 * plasma with iPIC3D." Mathematics and Computers in Simulation 80.7 (2010): 1509-1519.'
 *
 *        Copyright 2015 KTH Royal Institute of Technology
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at 
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mpi.h"
#include "ipicdefs.h"
#include "Basic.h"
#include "Alloc.h"
#include "TimeTasks.h"
#include "errors.h"

/** method to calculate extrapolation of 6 neighbours cells values to the outermost cell*/
double extrapolate(double f[7], double weig[3], int dir[6]) {
  double dxf, dyf, dzf,den,result;
  double f_closer=0.0; 
  int ndir[6];
  for (int i=0; i<6; i++){
    ndir[i]=1-dir[i];
    //f_closer+=(dir[i]*f[i+1]);
    f_closer+=f[i+1];
  }
  f_closer = f_closer/6.;
  //dxf = (f[1]*ndir[0]-f[2]*ndir[1])*weig[0];
  //dyf = (f[3]*ndir[2]-f[4]*ndir[3])*weig[1];
  //dzf = (f[5]*ndir[4]-f[6]*ndir[5])*weig[2];
  //den = 1.0+weig[0]*dir[0]-weig[0]*dir[1]+weig[1]*dir[2]-weig[1]*dir[3]+weig[2]*dir[4]-weig[2]*dir[5];
  //result = (f[0]+dxf+dyf+dxf)/den;
  result = f_closer;
  return result;
}

/** method to calculate the parallel dot product with vect1, vect2 having the ghost cells*/
double dotP(const double *vect1, const double *vect2, int n,MPI_Comm* comm) {
  double result = 0;
  double local_result = 0;
  for (register int i = 0; i < n; i++)
    local_result += vect1[i] * vect2[i];
  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, *comm);
  return (result);

}
/** method to calculate dot product */
double dot(const double *vect1, const double *vect2, int n) {
  double result = 0;
  for (int i = 0; i < n; i++)
    result += vect1[i] * vect2[i];
  return (result);
}
/** method to calculate the square norm of a vector */
double norm2(const double *const*vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      result += vect[i][j] * vect[i][j];
  return (result);
}
/** method to calculate the square norm of a vector */
double norm2(const arr3_double vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      result += vect.get(i,j,0) * vect.get(i,j,0);
  return (result);
}
/** method to calculate the square norm of a vector */
double norm2(const double *vect, int nx) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    result += vect[i] * vect[i];
  return (result);
}



/** method to calculate the parallel dot product */
/*
double norm2P(const arr3_double vect, int nx, int ny, int nz) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        local_result += vect.get(i,j,k) * vect.get(i,j,k);

  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (result);
}*/

/** method to calculate the parallel norm of a vector on different processors with the ghost cell */
/*
double norm2P(const double *vect, int n) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < n; i++)
    local_result += vect[i] * vect[i];
  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (result);
}*/
/** method to calculate the parallel norm of a vector on different processors with the gost cell*/
double normP(const double *vect, int n,MPI_Comm* comm) {
  double result = 0.0;
  double local_result = 0.0;
  for (register int i = 0; i < n; i++)
    local_result += vect[i] * vect[i];
  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, *comm);
  return (sqrt(result));

}
/** method to calculate the difference of two vectors*/
void sub(double *res, const double *vect1, const double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    res[i] = vect1[i] - vect2[i];
}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(double *vect1, const double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] += vect2[i];


}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) += vect2.get(i,j,k);
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr3_double vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) += vect2.get(i,j,0);
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr4_double vect2, int nx, int ny, int nz, int ns) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) += vect2.get(ns,i,j,k);
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr4_double vect2, int nx, int ny, int ns) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) += vect2.get(ns,i,j,0);
}
/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
void sub(arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) -= vect2.get(i,j,k);
}

/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
void sub(arr3_double vect1, const arr3_double vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) -= vect2.get(i,j,0);
}


/** method to sum 4 vectors vector1 = alfa*vector1 + beta*vector2 + gamma*vector3 + delta*vector4 */
void sum4(arr3_double vect1, double alfa, const arr3_double vect2, double beta, const arr3_double vect3, double gamma, const arr3_double vect4, double delta, const arr3_double vect5, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * (vect2.get(i,j,k) + beta * vect3.get(i,j,k) + gamma * vect4.get(i,j,k) + delta * vect5.get(i,j,k));

}
/** method to calculate the scalar-vector product */
void scale(double *vect, double alfa, int n) {
  for (register int i = 0; i < n; i++)
    vect[i] *= alfa;
}

/** method to calculate the scalar-vector product */
void scale(arr3_double vect, double alfa, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect.fetch(i,j,0) *= alfa;
}


/** method to calculate the scalar-vector product */
void scale(arr3_double vect, double alfa, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect.fetch(i,j,k) *= alfa;
}
/** method to calculate the scalar-vector product */
void scale(arr3_double vect1, const arr3_double vect2, double alfa, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k) * alfa;
}

/** method to calculate the scalar-vector product */
void scale(arr3_double vect1, const arr3_double vect2, double alfa, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) = vect2.get(i,j,0) * alfa;
}

/** method to calculate the scalar-vector product */
void scale(double *vect1, const double *vect2, double alfa, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] = vect2[i] * alfa;
}

/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(double alfa, arr3_double vect1, arr3_double vect2, const arr3_double vect3, int nx, int ny, int nz){
	  for (int i = 0; i < nx; i++)
	    for (int j = 0; j < ny; j++)
	      for (int k = 0; k < nz; k++)
	        vect3.fetch(i,j,k) = vect1.get(i,j,k) + alfa * vect2.get(i,j,k);
}
void addscale(double alfa, arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect1.get(i,j,k) + alfa * vect2.get(i,j,k);
}
/** add scale for weights */
void addscale(double alfa, double vect1[][2][2], double vect2[][2][2], int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect1[i][j][k] + alfa * vect2[i][j][k];

}
/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(double alfa, arr3_double vect1, const arr3_double vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) += alfa * vect2.get(i,j,0);
}
/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(double alfa, double *vect1, const double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] += alfa * vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2   */
void addscale(double alfa, double beta, double *vect1, const double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] = vect1[i] * beta + alfa * vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
void addscale(double alfa, double beta, arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz) {

  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++) {
        vect1.fetch(i,j,k) = beta * vect1.get(i,j,k) + alfa * vect2.get(i,j,k);
      }

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
void addscale(double alfa, double beta, arr3_double vect1, const arr3_double vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) = beta * vect1.get(i,j,0) + alfa * vect2.get(i,j,0);

}

/** method to calculate v1 = v1 + alpha * v2 * v3 */
void sumscalprod(double ***vect1, double alfa, double ***vect2, double ***vect3, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1[i][j][k] += alfa * vect2[i][j][k] * vect3[i][j][k];
}

/** method to calculate vector1 = alfa*vector2 + beta*vector3 */
void scaleandsum(arr3_double vect1, double alfa, double beta, const arr3_double vect2, const arr3_double vect3, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * vect2.get(i,j,k) + beta * vect3.get(i,j,k);
}
/** method to calculate vector1 = alfa*vector2 + beta*vector3 with vector2 depending on species*/
void scaleandsum(arr3_double vect1, double alfa, double beta, const arr4_double vect2, const arr3_double vect3, int ns, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * vect2.get(ns,i,j,k) + beta * vect3.get(i,j,k);
}
/** method to calculate vector1 = alfa*vector2*vector3 with vector2 depending on species*/
void prod(arr3_double vect1, double alfa, const arr4_double vect2, int ns, const arr3_double vect3, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * vect2.get(ns,i,j,k) * vect3.get(i,j,k);

}
/** method to calculate vect1 = vect2/alfa */
void div(arr3_double vect1, double alfa, const arr3_double vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k) / alfa;

}
void prod6(arr3_double vect1, const arr3_double vect2, const arr3_double vect3, const arr3_double vect4, const arr3_double vect5, const arr3_double vect6, const arr3_double vect7, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k) * vect3.get(i,j,k) + vect4.get(i,j,k) * vect5.get(i,j,k) + vect6.get(i,j,k) * vect7.get(i,j,k);
}
/** method used for calculating PI */
void proddiv(arr3_double vect1, const arr3_double vect2, double alfa, const arr3_double vect3, const arr3_double vect4, const arr3_double vect5, const arr3_double vect6, double beta, const arr3_double vect7, const arr3_double vect8, double gamma, const arr3_double vect9, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = (vect2.get(i,j,k) + alfa * (vect3.get(i,j,k) * vect4.get(i,j,k) - vect5.get(i,j,k) * vect6.get(i,j,k)) + beta * vect7.get(i,j,k) * vect8.get(i,j,k)) / (1 + gamma * vect9.get(i,j,k));

  // questo mi convince veramente poco!!!!!!!!!!!!!! CAZZO!!!!!!!!!!!!!!!!!!
  // ***vect1++ = (***vect2++ + alfa*((***vect3++)*(***vect4++) - (***vect5++)*(***vect6++)) + beta*(***vect7++)*(***vect8++))/(1+gamma*(***vect9++));
}
/** method to calculate the opposite of a vector */
void neg(arr3_double vect, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect.fetch(i,j,k) = -vect.get(i,j,k);
}

/** method to calculate the opposite of a vector */
void neg(arr3_double vect, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect.fetch(i,j,0) = -vect.get(i,j,0);
}
/** method to calculate the opposite of a vector */
void neg(arr3_double vect, int nx) {
  for (register int i = 0; i < nx; i++)
    vect.fetch(i,0,0) = -vect.get(i,0,0);
}
/** method to calculate the opposite of a vector */
void neg(double *vect, int n) {
  for (register int i = 0; i < n; i++)
    vect[i] = -vect[i];


}
/** method to set equal two vectors */
void eq(arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k);

}
/** method to set equal two vectors */
void eq(arr3_double vect1, const arr3_double vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) = vect2.get(i,j,0);

}

/** method to set equal two vectors */
void eq(arr4_double vect1, const arr3_double vect2, int nx, int ny, int is) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(is,i,j,0) = vect2.get(i,j,0);

}
/** method to set equal two vectors */
void eq(arr4_double vect1, const arr3_double vect2, int nx, int ny, int nz, int is) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(is,i,j,k) = vect2.get(i,j,k);

}

/** method to set a vector to a Value */
void eqValue(double value, arr3_double vect, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect.fetch(i,j,k) = value;

}
//void eqValue(double value, double vect[][2][2], int nx, int ny, int nz) {
//  for (int i = 0; i < nx; i++)
//    for (int j = 0; j < ny; j++)
//      for (int k = 0; k < nz; k++)
//        vect[i][j][k] = value;
//
//}
/** method to set a vector to a Value */
void eqValue(double value, arr3_double vect, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect.fetch(i,j,0) = value;

}
/** method to set a vector to a Value */
void eqValue(double value, arr3_double vect, int nx) {
  for (register int i = 0; i < nx; i++)
    vect.fetch(i,0,0) = value;

}
/** method to set a vector to a Value */
void eqValue(double value, double *vect, int n) {
  for (register int i = 0; i < n; i++)
    vect[i] = value;
}
/** method to put a column in a matrix 2D */
void putColumn(double **Matrix, double *vect, int column, int n) {
  for (int i = 0; i < n; i++)
    Matrix[i][column] = vect[i];

}
/** method to get a column in a matrix 2D */
void getColumn(double *vect, double **Matrix, int column, int n) {
  for (int i = 0; i < n; i++)
    vect[i] = Matrix[i][column];
}
/** method to get rid of the ghost cells */
void getRidGhost(double **out, double **in, int nx, int ny) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      out[i - 1][j - 1] = in[i][j];
}
