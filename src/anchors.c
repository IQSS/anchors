/*

  Functions: vecrowsum, opll, opllgh
  Author   : Jonathan Wand <jwand@latte.harvard.edu>

  Purpose  : evaluation of ordered probit, ordered probit with RE, and summing rows of a matrix

  Created:   2002-12-27
  Modified:  $Date: 2005/08/10 23:25:49 $
  Revision:  $Revision: 1.3 $
  RCS-ID:    $Id: anchors.c,v 1.3 2005/08/10 23:25:49 jwand Exp $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

extern double pnorm(double x,double mean,double sd,int lower_tail,int log_p);
extern double plnorm(double x,double logmean,double logsd,int lower_tail,int log_p);
extern double dnorm(double x,double mean,double sd,int log_p);



/* 
   For both opll and opllgh 

   y    : real vector of choices (1..ncats), self1 stacked above self2 above self3 etc 
   taus : real vector of cumulative taus, self1-cut1 above self1-cut2... self2-cut1 above self2-cut2 
   xb   : real vector of means 

*/

void opllgrself( int *xnobs, 
		 int *xncat,  
		 int *xnvarx,  
		 int *xnvarv,  
		 int *xnvarv1,  
		 int *xnself,
		 double *xse, 
		 int   *y, 
		 double *xb, 
		 double *taus, 
		 double *V,
		 double *V1,
		 double *X,
		 double *dLdSigma, /* scaler */ 
		 double *dLdBeta,   /* 1 x nvar.x */
		 double *dLdGamma,  /* 1 x nvar.v */
		 double *dLdGamma1  /* 1 x nvar.v */
		 ) ;

void opllgrvign( int *xnobs, 
		 int *xncat,  
		 //int *xnvarx,  
		 int *xnvarv,  
		 int *xnvarv1,  
		 double *xse, 
		 int   *y, 
		 double *xb, 
		 double *taus, 
		 double *V,
		 double *V1,
		 //double *X,
		 double *dLdSigma, /* scaler */ 
		 double *dLdBeta,   /* 1 x nvar.x */
		 double *dLdGamma,  /* 1 x nvar.v */
		 double *dLdGamma1  /* 1 x nvar.v */
		 ) ;

void opll( int *xnobs, int *xncat,  
	   double *xse, double *xpenalty, 
	   int *y, 
	   double *xb, double *taus, 
	   double *llik);

void opllgh( int *xnobs, int *xngh, int *xnself, int *xncat,  
	     double *xs2su, double *xse, double *xpenalty, 
	     double *ghx, double *ghw, 
	     int *y, 
	     double *xb, double *taus, 
	     double *llik);

void opllgh( int *xnobs, int *xngh, int *xnself, int *xncat,  
	     double *xs2su, double *xse, double *xpenalty, 
	     double *ghx, double *ghw, 
	     int *y, 
	     double *xb, double *taus, 
	     double *llik) {

  int nobs=*xnobs, ngh=*xngh, nself=*xnself, ncat=*xncat;
  double s2su = *xs2su, se = *xse, penalty = *xpenalty;

  int i,j,m;
  int choice;
  double p, u;
  double *mfjm ; 

  mfjm = (double *) calloc( (int) nobs,sizeof(double));
  for (i=0; i<nobs; i++) {
    mfjm[i] = 0.0;
  }
  /*
  printf("In opllgh\n");
  fflush(stdout);
  */
  /*
  printf("nobs %d ngh %d nself %d ncat %d\n",
	 nobs,ngh,nself,ncat);
  printf("s2su %g se %g penalty %g\n",
	 s2su,se,penalty);
  fflush(stdout);
  */

  /* 
     Integrate over random effect, by evaluating at each m nodes of GH
     Evaluate only choices actually used by an individual
  */ 
  for (m=0; m < ngh; m++) {
    u = s2su * ghx[m];

    /* LOOP over individuals 'i' */
    for (i=0; i<nobs; i++) {

      p = 1.0;
      for (j=0; j < nself; j++) {
	choice = y[ j*nobs + i ];
	/*
	printf("m %d i %d j %d choice %d index-1 %d index-2 %d\n",
	       m,i,j,choice, 
	       (j * nobs * (ncat-1)) + (choice-1)*nobs + i,
	       (j * nobs * (ncat-1)) + (choice-2)*nobs + i);
	fflush(stdout);  
	*/
	if (choice==1) {
	  /* 
	  printf("p1 %g\n",
		 taus[ (j * nobs * ncat) + (choice-1)*nobs + i ] - xb[i]);
	  fflush(stdout);
	  */
	  p *= pnorm( (taus[ (j * nobs * (ncat-1)) + (choice-1)*nobs + i ] - xb[i]) / se - u ,0.0,1.0,1,0);
	}
	else if ((choice > 1) & (choice < ncat)) {
	  /*
	  printf("p2 taus %g - %g\n",
		 taus[ (j * nobs * ncat) + (choice-1)*nobs + i ] - xb[i],
		 taus[ (j * nobs * ncat) + (choice-2)*nobs + i ] - xb[i]);
	  printf("se %g u %g\n",se,u);
	  fflush(stdout); 
	  */
	  p *= pnorm( (taus[ (j * nobs * (ncat-1)) + (choice-1)*nobs + i ] - xb[i]) / se - u ,0.0,1.0,1,0) -
                pnorm( (taus[ (j * nobs * (ncat-1)) + (choice-2)*nobs + i ] - xb[i]) / se - u ,0.0,1.0,1,0);
	}
	else if ( choice == ncat ) {
	  /*
	  printf("p3 %g\n",
		 taus[ (j * nobs * ncat) + (choice-2)*nobs + i ] - xb[i]);
	  fflush(stdout); 
	  */
	  p *= pnorm((taus[ (j * nobs * (ncat-1)) + (choice-2)*nobs + i ] - xb[i]) / se - u ,0.0,1.0,0,0);
	}
	else {
	  /*
	  printf("p4\n");
	  fflush(stdout);
	  */
	  /* do nothing -- does not contribute to LL */
	  /* if all self-questions for person 'i' are missing, then p=1 */
	}
      }
      /*
      printf("prob %g\n",p);
      fflush(stdout); 
      */
      mfjm[i] +=  ghw[m] * p;    
      /*
      printf("mjfm %g\n",mfjm[i]);
      fflush(stdout);
      */
    }
  }
  /*
  printf("get llik\n");
  fflush(stdout);
  */

  /* do error checking... penalize with 'penalty' non-valid observations */
  /* NOTE M_SQRT_PI = 1.772453850905516027298167483341 */
  for (i=0; i < nobs; i++) {
    /*
    printf("ii %d\n",i);
    fflush(stdout);
    */
    if ( ISNA(mfjm[i]) ||  mfjm[i] <= 0 || mfjm[i] > 1) 
      llik[i] = penalty;
    else
      llik[i] = log( mfjm[i] / M_SQRT_PI );
    /*
    printf("i %d llik %g mfjm %g penalty %g sqrtpi %g\n",
	   i, llik[i], mfjm[i],
	   penalty,M_SQRT_PI);
    */
  }




  free(mfjm);

  /*
  printf("Leaving opllgh\n");
  fflush(stdout);
  */
}



void opll( int *xnobs, int *xncat,  
	     double *xse, double *xpenalty, 
	     int *y, 
	     double *xb, double *taus, 
	     double *llik) {

  int nobs=*xnobs, ncat=*xncat;
  double se = *xse, penalty = *xpenalty;

  int i;
  int choice;


  /* LOOP over individuals 'i' */
  for (i=0; i<nobs; i++) {
    choice = y[ i ];

    if (choice==1) {

      llik[i] = pnorm( (taus[ (choice-1)*nobs + i ] - xb[i]) / se ,0.0,1.0,1,0);
    }
    else if ((choice > 1) & (choice < ncat)) {
      llik[i] = pnorm( (taus[ (choice-1)*nobs + i ] - xb[i]) / se ,0.0,1.0,1,0) -
                pnorm( (taus[ (choice-2)*nobs + i ] - xb[i]) / se ,0.0,1.0,1,0);
    }
    else if ( choice == ncat ) {
      llik[i] = pnorm( (taus[ (choice-2)*nobs + i ] - xb[i]) / se ,0.0,1.0,0,0);
    }
    else {
      /* do nothing -- does not contribute to LL */
      /* if all self-questions for person 'i' are missing, then p=1 */
      llik[i] = 1;
    }
  }

  /* do error checking... penalize with 'penalty' non-valid observations */
  for (i=0; i < nobs; i++) {
    /*
    printf("ii %d\n",i);
    fflush(stdout);
    */
    if ( ISNA(llik[i]) ||  llik[i] <= 0 || llik[i] > 1) 
      llik[i] = penalty;
    else
      llik[i] = log(llik[i]);
  }
}

void opllgrself( int *xnobs, 
		 int *xncat,  
		 int *xnvarx,  
		 int *xnvarv,  
		 int *xnvarv1,  
		 int *xnself,
		 double *xse, 
		 int   *y, 
		 double *xb, 
		 double *taus, 
		 double *V,
		 double *V1,
		 double *X,
		 double *dLdSigma, /* scaler */ 
		 double *dLdBeta,  /* 1 x nvar.x */
		 double *dLdGamma,  /* 1 x nvar.v * n.cat */
		 double *dLdGamma1  /* 1 x nvar.v * n.cat */
		 
		 ) {

  int nobs=*xnobs, ncat=*xncat, nvarx = *xnvarx, nvarv = *xnvarv, nvarv1 = *xnvarv1,
       nself = *xnself;  
  double se = *xse;

  int i,choice,k,m,idx2,idxK, idx1,idxM,self; 
  double tmp;
  double m1, m2, denom;

  double pmat1,pmat2,dmat1,dmat2,vmat1,vmat2,tmp1,tmp2;

  for (i=0; i < nobs; i++) {
    
    for (self=0; self<nself; self++) {
      choice  = y[self*nobs+i]; /* choices indexed 1:ncat */
      if (choice > 0) {

        if (choice==1) {
	  idx1 = -99;
	  tmp1 = -99;
	  idx2 = (self * nobs * (ncat-1)) + (choice-1)*nobs + i;
	  tmp2 = (taus[ idx2  ] - xb[i]) / se ;
	  vmat1 = 0;
	  dmat1 = 0;
	  pmat1 = 0;
	  vmat2 = tmp2;
	  dmat2 = dnorm( tmp2 ,0.0,1.0,0);
	  pmat2 = pnorm( tmp2 ,0.0,1.0,1,0);
        } 
        else if (choice==ncat) {
	  idx1 = (self * nobs * (ncat-1)) + (choice-2)*nobs + i;
	  tmp1 = (taus[ idx1 ] - xb[i]) / se ;
	  idx2 = -99;
	  tmp2 = -99;
	  vmat1 = tmp1;
	  dmat1 = dnorm( tmp1 ,0.0,1.0,0);
	  pmat1 = pnorm( tmp1 ,0.0,1.0,1,0);
	  vmat2 = 0;
	  dmat2 = 0;
	  pmat2 = 1;
        }
        else {
	  idx1 = (self * nobs * (ncat-1)) + (choice-2)*nobs + i;
	  tmp1 = (taus[ idx1 ] - xb[i]) / se ;
	  idx2 = (self * nobs * (ncat-1)) + (choice-1)*nobs + i;
	  tmp2 = (taus[ idx2  ] - xb[i]) / se ;
	  vmat1 = tmp1;
	  dmat1 = dnorm( tmp1 ,0.0,1.0,0);
	  pmat1 = pnorm( tmp1 ,0.0,1.0,1,0);
	  vmat2 = tmp2;
	  dmat2 = dnorm( tmp2 ,0.0,1.0,0);
	  pmat2 = pnorm( tmp2 ,0.0,1.0,1,0);
        }
//	printf("i %d %d : %g %g : %g %g : %g %g : %g %g \n",
//	       i,choice,
//	       pmat1,pmat2,dmat1,dmat2,vmat1,vmat2,tmp1,tmp2);
//	fflush(stdout);
//	printf("i %d %d : %d %d : %g %g : %g %g : %g %g : %g %g : xb %g \n",
//	       i,choice,idx1,idx2,
//	       pmat1,pmat2,vmat1,vmat2,tmp1,tmp2,vmat1+xb[i],vmat2+xb[i],xb[i]);
//	fflush(stdout);


	denom = pmat2 - pmat1;
//	printf("denom %g \n",denom);
//	fflush(stdout);
  
        tmp = (dmat1 * vmat1 - dmat2 * vmat2)/ denom;
        *dLdSigma = *dLdSigma + tmp;
//	printf("S %g\n",*dLdSigma);
//	fflush(stdout);
    
        for (k=0; k<nvarx; k++) {
          idxK = k * nobs + i;
          tmp  = ( dmat1 - dmat2) / denom;
          dLdBeta[k] = dLdBeta[k] + (X[idxK] * tmp);
        }
    
        m1 = (choice == ncat ? 0.0 : 1.0);


        for (k=0; k<1; k++) {
          idxK = self*(ncat-1)*nvarv1 + k*nvarv1; /* still need to add 'm' to this */
    
          for (m=0; m<nvarv1; m++) {
	    idxM =  m * nobs + i;
    
    
	    m2 = (k <= (choice-2) ? 1.0 : 0.0); 
    
	    if (k <= (choice-1)) {
	      dLdGamma1[idxK+m] = 
		dLdGamma1[idxK+m] +
		(m1 * V1[idxM] * dmat2 -
		 m2 * V1[idxM] * dmat1) / denom;

	    }
          }
	}

        for (k=1; k<(ncat-1); k++) {
          idxK = self*(ncat-1)*nvarv + (k-1)*nvarv; /* still need to add 'm' to this */
    
          for (m=0; m<nvarv; m++) {
	    idxM =  m * nobs + i;
    
    
	    m2 = (k <= (choice-2) ? 1.0 : 0.0); 
    
	    if (k <= (choice-1)) {
	      dLdGamma[idxK+m] = 
		dLdGamma[idxK+m] +
		(m1 * V[idxM] * dmat2 -
		 m2 * V[idxM] * dmat1) / denom;

//	      printf("i %d c %d k %d m %d : idxK+m %d %d : dG %g : V %g\n",
//		     i,choice,k,m,idxK+m,idxM,dLdGamma[idxK+m],V[idxM]);

	    }
    
          }
        }
      }
    }
  }

  *dLdSigma = *dLdSigma / se;
  for (k=0; k<nvarx; k++) {
    dLdBeta[k] =  dLdBeta[k] / se;
  }
  for (self=0; self<nself; self++) {
    for (k=0; k<1; k++) {
      idxK = self*(ncat-1)*nvarv1 + k*nvarv1; 
      for (m=0; m<nvarv1; m++) {
	dLdGamma1[idxK+m] = dLdGamma1[idxK+m] / se;
      }
    }
    for (k=1; k<(ncat-1); k++) {
      idxK = self*(ncat-1)*nvarv + (k-1)*nvarv; 
      for (m=0; m<nvarv; m++) {
	dLdGamma[idxK+m] = dLdGamma[idxK+m] / se;
      }
    }
  }

}

void opllgrvign( int *xnobs, 
		 int *xncat,  
		 //int *xnvarx,  
		 int *xnvarv,  
		 int *xnvarv1,  
		 double *xse, 
		 int   *y, 
		 double *xb, 
		 double *taus, 
		 double *V,
		 double *V1,
		 //double *X,
		 double *dLdSigma, /* scaler */ 
		 double *dLdTheta,  /* 1 x nvar.x */
		 double *dLdGamma,  /* 1 x nvar.v * n.cat */
		 double *dLdGamma1  /* 1 x nvar.v * n.cat */
		 
		 ) {

  int nobs=*xnobs, ncat=*xncat, nvarv = *xnvarv, nvarv1 = *xnvarv1;  
  double se = *xse;

  int i,choice,k,m,idx2,idxK, idx1,idxM; 
  double tmp;
  double m1, m2, denom;

  double pmat1,pmat2,dmat1,dmat2,vmat1,vmat2,tmp1,tmp2;

  for (i=0; i < nobs; i++) { // loop over individuals

    choice  = y[i]; // choices indexed 1:ncat 

    if (choice > 0) { // only valid cases


//      printf("iv %d choice %d : %d %d : %g %g : %g %g : %g %g \n",
//	     i,choice,idx1,idx2,tmp1,tmp2,taus[idx1],taus[idx2],xb[i],se);
//      fflush(stdout);

      if (choice==1) {
	idx2  =  (choice-1) * nobs + i   ;
	tmp2 = (taus[ idx2  ] - xb[i]) / se ;
	vmat1 = 0;
	dmat1 = 0;
	pmat1 = 0;
	vmat2 = tmp2;
	dmat2 = dnorm( tmp2 ,0.0,1.0,0);
	pmat2 = pnorm( tmp2 ,0.0,1.0,1,0);
      } 
      else if (choice==ncat) {
	idx1 =  (choice-2) * nobs + i   ;
	tmp1 = (taus[ idx1 ] - xb[i]) / se ;
	vmat1 = tmp1;
	dmat1 = dnorm( tmp1 ,0.0,1.0,0);
	pmat1 = pnorm( tmp1 ,0.0,1.0,1,0);
	vmat2 = 0;
	dmat2 = 0;
	pmat2 = 1;
      }
      else {
	idx1 =  (choice-2) * nobs + i   ;
	tmp1 = (taus[ idx1 ] - xb[i]) / se ;
	idx2  =  (choice-1) * nobs + i   ;
	tmp2 = (taus[ idx2  ] - xb[i]) / se ;
	vmat1 = tmp1;
	dmat1 = dnorm( tmp1 ,0.0,1.0,0);
	pmat1 = pnorm( tmp1 ,0.0,1.0,1,0);
	vmat2 = tmp2;
	dmat2 = dnorm( tmp2 ,0.0,1.0,0);
	pmat2 = pnorm( tmp2 ,0.0,1.0,1,0);
      }
      //      printf("%g %g : %g %g : %g %g\n",pmat1,pmat2,dmat1,dmat2,vmat1,vmat2);
      //      fflush(stdout);


      denom = pmat2 - pmat1;
      //      printf("denom %g \n",denom);
      //      fflush(stdout);

      tmp = (dmat1 * vmat1 - dmat2 * vmat2)/ denom;
      *dLdSigma = *dLdSigma + tmp;
  

      tmp  = ( dmat1 - dmat2) / denom;
      *dLdTheta = *dLdTheta + tmp;
  
      m1 = (choice == ncat ? 0.0 : 1.0);
      for (k=0; k<1; k++) {
        idxK = k*nvarv1; /* still need to add 'm' to this */
  
	//	printf("G %d ",k);
        for (m=0; m<nvarv1; m++) {
	  idxM =  m * nobs + i;
  
  
	  m2 = (k <= (choice-2) ? 1.0 : 0.0); 
  
	  if (k <= (choice-1)) {
	    dLdGamma1[idxK+m] = dLdGamma1[idxK+m] + 
	      (m1 * V1[idxM] * dmat2 -
	       m2 * V1[idxM] * dmat1) / denom;
	    //	    printf("%g ",dLdGamma[idxK+m]);
	  }
        }
	//printf("\n");
      }

      for (k=1; k<(ncat-1); k++) {
        idxK = (k-1)*nvarv; /* still need to add 'm' to this */
  
	//	printf("G %d ",k);
        for (m=0; m<nvarv; m++) {
	  idxM =  m * nobs + i;
  
  
	  m2 = (k <= (choice-2) ? 1.0 : 0.0); 
  
	  if (k <= (choice-1)) {
	    dLdGamma[idxK+m] = dLdGamma[idxK+m] + 
	      (m1 * V[idxM] * dmat2 -
	       m2 * V[idxM] * dmat1) / denom;
	    //	    printf("%g ",dLdGamma[idxK+m]);
	  }
        }
	//printf("\n");
      }
//      printf("S %g\n",*dLdSigma);
//      printf("T %g\n",*dLdTheta);
//      fflush(stdout);
    }
  }

  *dLdSigma = *dLdSigma / se;
  *dLdTheta = *dLdTheta / se;

  for (k=0; k<1; k++) {
    idxK = k*nvarv1; 
    for (m=0; m<nvarv1; m++) {
      idxM =  m * nobs + i;
      dLdGamma1[idxK+m] = dLdGamma1[idxK+m] / se;
    }
  }
  for (k=1; k<(ncat-1); k++) {
    idxK = (k-1)*nvarv; 
    for (m=0; m<nvarv; m++) {
      idxM =  m * nobs + i;
      dLdGamma[idxK+m] = dLdGamma[idxK+m] / se;
    }
  }
  //  printf("EXIT\n");
}


