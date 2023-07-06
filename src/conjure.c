/*

  Functions: vecrowsum
  Author   : Jonathan Wand <jwand@latte.harvard.edu>

  Purpose  : summing rows of a matrix

  Created:   2002-12-27
  Modified:  $Date:  $
  Revision:  $Revision: $
  RCS-ID:    $Id: $

*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include "Rmath.h"
void vecrowcumsum( double *taus, int *xnself, int *xncat, int *xnobs) ;

void vecrowcumsum( double *taus, int *xnself, int *xncat, int *xnobs) {

  int i,j,k,offset;
  int ncat = *xncat, nself = *xnself, nobs = *xnobs; 

  /*
  printf("ncat %d nself %d nobs %d\n",ncat,nself,nobs);
  fflush(stdout);
  */

  /* LOOP over individuals 'i' */
  for (j=0; j < nself; j++) {
    offset = j * nobs * (ncat);
    for (k=1; k<(ncat); k++) {
      for (i=0; i<nobs; i++) {

	/*
	printf("offset %d k %d i %d index %d index2 %d \n", offset, k, i ,
	       offset + k * nobs + i, 
	       offset + (k-1) * nobs + i);
	fflush(stdout);
	*/

	taus[ offset + k * nobs + i ] += taus[ offset + (k-1)* nobs + i ] ;
      }
    }
  }
}

void diff_nonneg( double *y, int *xn, int *xk) ;

void diff_nonneg( double *y, int *xn, int *xk) {
  int i;
  int n = *xn, k = *xk;

  for (i=0; i< (n-k); i++) {
    y[i] = ( y[i+k] >= y[i] ? 1 : 0);    /* get T/F for positive difference */
  }
}

void diff_sign( double *y, int *xn, int *xk) ;

void diff_sign( double *y, int *xn, int *xk) {
  int i, tmp;
  int n = *xn, k = *xk;

  for (i=0; i< (n-k); i++) {
    tmp = y[i+k] - y[i];
    y[i] = ( tmp >= 0 ? 1 : 0);    /* get signed of difference */
    y[i] = ( tmp == 0 ? 0 : y[i]); /* ties are equal to zero   */
  }
}

void Tcount( double *pos, int *xk, int *xn, int *xrA, int *xrB, double *out) ;

void Tcount( double *pos, int *xk, int *xn, int *xrA, int *xrB, double *out) {
  int r,i,j;
  int offset = 0;
  int n = *xn, k = *xk, rA = *xrA, rB = *xrB ;
  double sum, weight, half;
  
  /* loop over subinterval widths */
  for (r = rA; r <= rB; r++) {
    weight = 1/sqrt(r);
    half   =  r/2.0;
    /* loop over starting points */
    for (j = 0; j < (n-r-k+1) ; j++ ) {
      /* sum indicators for positive difference in subinterval*/
      sum = 0.0;
      for (i = j; i < j+r; i++) {
	sum = sum+pos[i];
//	printf("r=%d j=%d i=%d \n",
//	       r,j,i);
//	fflush(stdout);
      }
      /* and record weighted difference from half of r */
      out[ offset + j ]  = weight * ( half - sum );
//      printf("%d %d %d : %g = %g * (%g  - %g)\n",
//	     r,j,offset+j,
//	     out[ offset + j ],
//	     weight,half,sum);
//      fflush(stdout);

    } 
    offset += (n-r-k+1); /* increment offset for next round */
  }
}

