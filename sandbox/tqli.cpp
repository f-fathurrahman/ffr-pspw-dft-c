#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SQR(x) ((x)*(x))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define IDX2F(i,j,DIM1) (((j)-1)*(DIM1) + ((i)-1))

double pythag(double a, double b)
// Computes (a**2 + b**2)**(1/2) without destructive underflow or overflow.
{
  double absa,absb;
  absa = fabs(a);
  absb = fabs(b);
  if(absa > absb) return absa*sqrt(1.0 + SQR(absb/absa));
  else return(absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa/absb)));
}

void eigsrt(double *d, double *v, int n)
/* Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi
 * tqli, this routine sorts the eigenvalues into descending order, and rearranges
 * the columns of v correspondingly. The method is straight insertion.
 */
{
  int k,j,i;
  double p;
  for(i=1; i<n; i++) {
    k = i;
    p = d[k-1];
    for (j=i+1; j<=n; j++) {
      if (d[j-1] >= p) {
        k = j;
        p = d[k-1];
      }
    }
    if (k != i) {
      d[k-1] = d[i-1];
      d[i-1] = p;
      for (j=1; j<=n; j++) {
        p = v[IDX2F(j,i,n)];
        v[IDX2F(j,i,n)] = v[IDX2F(j,k,n)];
        v[IDX2F(j,k,n)] = p;
      }
    }
  }
}

void tqli(double *d, double *e, int n, double *z)
/* QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
 * symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2
 * x11.2. On input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On
 * output, it returns the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of
 * the tridiagonal matrix, with e[1] arbitrary. On output e is destroyed. When finding only the
 * eigenvalues, several lines may be omitted, as noted in the comments. If the eigenvectors of
 * a tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as the identity matrix.
 * If the eigenvectors of a matrix that has been reduced by tred2 are required, then z is input
 * as the matrix output by tred2. 
 * In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].
 */
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  for(i=2; i<=n; i++) e[i-2] = e[i-1]; // Convenient to renumber the elements of e. e[n]=0.0;

  for(l=1; l<=n; l++) {
    iter=0;
    do {
      // Look for a single small subdiagonal element to splitthe matrix.
      for (m=l; m<=n-1; m++) { 
        dd = fabs(d[m-1]) + fabs(d[m]);
        if ((double) (fabs(e[m-1]) + dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("ERROR: Too many iterations in tqli\n");
          abort();
        }
        g = (d[l] - d[l-1])/(2.0*e[l-1]); // Form shift.
        r = pythag(g,1.0);
        g = d[m-1] - d[l-1] + e[l-1]/(g + SIGN(r,g)); // This is dm - ks.
        s = c = 1.0;
        p = 0.0;
// A plane rotation as in the original QL, followed by Givens rotations
// to restore tridiagonal form.
        for(i=m-1; i>=l; i--) { 
          f = s*e[i-1];
          b = c*e[i-1];
          e[i] = ( r = pythag(f,g) );
          // Recover from underflow
          if (r == 0.0) {
            d[i] -= p;
            e[m-1] = 0.0;
            break;
          }
          s = f/r;
          c = g/r;
          g = d[i] - p;
          r = (d[i-1]-g)*s + 2.0*c*b;
          d[i] = g + (p=s*r);
          g = c*r-b;
          /* Next loop can be omitted if eigenvectors not wanted*/
          for (k=1;k<=n;k++) { // Form eigenvectors.
            f = z[IDX2F(k,i+1,n)];
            z[IDX2F(k,i+1,n)] = s*z[IDX2F(k,i,n)] + c*f;
            z[IDX2F(k,i,n)] = c*z[IDX2F(k,i,n)] - s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l-1] -= p;
        e[l-1] = g;
        e[m-1] = 0.0;
      }
    } while (m != l);
  }
}

