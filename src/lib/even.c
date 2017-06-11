/*
      harmonic1.c
  
        Solution of the quantum harmonic oscillator
        Eigenvalue search using shooting method.
        Forward and backward integration with Numerov method.
        Solution matching at a classical turning point.
        Adimensional units: x = (mK/hbar^2)^(1/4) X
                            e = E/(hbar omega)
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "config.h"

static double pi = 3.14159265358979323846;

typedef struct Outfile {
	FILE *y;
//	FILE *x, *p, *vpot, *f, *y_sq;
	FILE *e;
} Outfile;

int closeAllOutFile(Outfile out) {
	fclose(out.y);
//	fclose(out.x); fclose(out.p);
//	fclose(out.vpot); fclose(out.f); fclose(out.y_sq);
	fclose(out.e);
	return 0;
}

int main()
{
    double sqrt();

    int i, icl = -1; // mesh
    int hnodes = -1, ncross, parity, kkk, n_iter; // nodes
    double dx, ddx12, xmcl, norm, arg, yicl, djump; // xmax
    double elw, eup, e;
    double *x, *y, *p, *vpot, *f;
//  char fileout[80];
//	FILE *out;

	Outfile OutFiles; // Declaring bunch of output files

/* 	Process Config file */
	Config_t Config;
	processConfigFile(&Config, "./CONFIG");

    e = Config.e;

/*  Read input data */

//  fprintf(stderr, "Max value of x (typical value: 10) ? ");
//  scanf("%lf",&xmax);
//  fprintf(stderr, "Number of grid points (typically a few hundreds) ? " );
//  scanf("%d",&Config.mesh);

/*  Allocate arrays (from 0 to Config.mesh), Initialize grid */

    x = (double *) malloc( (Config.mesh+1) * sizeof (double));
    y = (double *) malloc( (Config.mesh+1) * sizeof (double));
    p = (double *) malloc( (Config.mesh+1) * sizeof (double));
    f = (double *) malloc( (Config.mesh+1) * sizeof (double));
    vpot = (double *) malloc( (Config.mesh+1) * sizeof (double));
    dx = Config.xmax / Config.mesh;
    ddx12 = dx * dx / 12.;


/*  set up the potential (must be even w.r.t. x=0) */

    for (i = 0; i <= Config.mesh; ++i) {
	x[i] = (double) i * dx;
	//vpot[i] = 0.5 * x[i] * x[i];
	/* Double well potential */
	vpot[i] = x[i]*x[i]*x[i]*x[i]*0.2*0.2*0.2*0.2 - 2*x[i]*x[i]*0.2*0.2 + 1;
    }
//  fprintf(stderr, "Output file name = ");
//  scanf("%80s", fileout);
//	out = fopen(Config.fileout, "w");

	//fprintf(stderr,"Output file writing begins.\n");
//	OutFiles.x = fopen("x.dat","w");
	OutFiles.y = fopen("y.dat","w");
//	OutFiles.y_sq = fopen("y_sq.dat","w");
//	OutFiles.p = fopen("p.dat","w");
//	OutFiles.vpot = fopen("vpot.dat","w");
//	OutFiles.f = fopen("f.dat","w");
	OutFiles.e = fopen("energy.dat","w");
	//fprintf(stderr,"Output file writing has been done.\n");

	/* writing x coord for each files */
	for (i = Config.mesh; i >= 1; --i) {
		fprintf(OutFiles.y, "%7.3f ",-x[i]);
	}
   	for (i = 0; i <= Config.mesh; ++i) {
		fprintf(OutFiles.y, "%7.3f ",x[i]);
   	}
	fprintf(OutFiles.y,"\n");

	
	int calculatedEigenVecNum = 0;

    /* this is the entry point for a new eigenvalue search */
    for (Config.nodes = Config.nodesStart; Config.nodes <= Config.nodesEnd; Config.nodes++) {
      
      if (!(calculatedEigenVecNum % 10)) fprintf(stderr,"Calculated EigenVecNum: %d\n",calculatedEigenVecNum);
	
	/* Initialization per eigenstate calc. */
	  if (Config.e == 0) {e = 0;}
	  
      /*  Read number of nodes (stop if < 0) */
//    fprintf(stderr, "Number of nodes (-1=exit) ? ");
//    scanf("%d",&nodes);

	  // ADDED .. caculatedEigenVecNum .. should be modified
	  /* 170519 modified */
	  /*
      if (Config.nodes < 0 || calculatedEigenVecNum >= 1) {
        free(vpot); free(f); free(p); free(y); free(x);
        fclose(out);
        exit(0);
      }*/
      
      /*  set initial lower and upper bounds to the eigenvalue */
      
      eup = vpot[Config.mesh];
      elw = eup;
      for (i = 0; i <= Config.mesh; ++i) {
        if ( vpot[i] < elw )
	  elw = vpot[i];
        if ( vpot[i] > eup )
	  eup = vpot[i];
      }
      
      /*  set trial energy */
      
//    fprintf(stderr, "Trial energy (0=search with bisection) ? ");
//    scanf("%lf", &e);
      if (e == 0.) { /* search eigenvalues with bisection (max 1000 iterations) */
	e = 0.5 * (elw + eup);
	n_iter = 1000;
      } else {	   /*  test a single energy value */
	n_iter = 1;
      }
      
      for (kkk = 0; (kkk < n_iter) && (eup-elw > 1.e-10); kkk++) {
	
	/*
	  set up the f-function used by the Numerov algorithm
	  and determine the position of its last crossing, i.e. change of sign
	  f < 0 means classically allowed   region
	  f > 0 means classically forbidden region
	*/
	f[0] = ddx12 * (2.*(vpot[0] - e));
	icl = -1;
	for (i = 1; i <= Config.mesh; ++i) {
          f[i] = ddx12 * 2. * (vpot[i] - e);
	  /*
	    beware: if f(i) is exactly zero the change of sign is not observed
	    the following line is a trick to prevent missing a change of sign 
	    in this unlikely but not impossible case:
	  */
	  if (f[i] == 0.) {
	    f[i] = 1e-20;
	  }
	  /*   store the index 'icl' where the last change of sign has been found */
	  if (f[i] != copysign(f[i],f[i-1])) {
	    icl = i;
	  }
	}
//	fprintf(stderr,"%d-th index: icl == %d\n",i,icl);	
	if (icl >= Config.mesh - 2) {
          fprintf(stderr, "last change of sign too far.");
          exit(1);
	}
	if (icl < 1) {
          fprintf(stderr, "no classical turning point?");
          exit(1);
	}
	
	/*   f(x) as required by the Numerov algorithm  */
	
	for (i = 0; i <= Config.mesh; ++i) {
          f[i] = 1. - f[i];
	}
	
	for (i = 0; i <= Config.mesh; ++i) {
          y[i] = 0.;
	}
	
	/*  Determination of the wave-function in the first two points  */
	
	hnodes = Config.nodes / 2;
	
	/*  beware the integer division: 1/2 = 0 !
	    if nodes is even, there are 2*hnodes nodes
	    if nodes is odd,  there are 2*hnodes+1 nodes (one is in x=0)
	    hnodes is thus the number of nodes in the x>0 semi-axis (x=0 excepted) */
	
	if (2*hnodes == Config.nodes) {
          /*  even number of nodes: wavefunction is even  */
	  y[0] = 1.;
          /*  assume f(-1) = f(1)  */
	  y[1] = 0.5 * (12. - f[0] * 10.) * y[0] / f[1];
	} else {
          /*  odd  number of nodes: wavefunction is odd   */
	  y[0] = 0.;
	  y[1] = dx;
	}
	
	/*   Outward integration and count number of crossings  */
	
	ncross = 0;
	for (i = 1; i <= icl-1; ++i) {
          y[i + 1] = ((12. - f[i] * 10.) * y[i] - f[i - 1] * y[i - 1])
	    / f[i + 1];
          if (y[i] != copysign(y[i],y[i+1]))
	    ++ncross;
	}
	yicl = y[icl];
	
	if (2*hnodes == Config.nodes) {
          /* even number of nodes: no node in x=a 0*/
          ncross = 2*ncross;
	} else {
          /*  odd number of nodes:    node in x=0 */
          ncross = 2*ncross+1;
	}
	
	/*  Check number of crossings  */
	
	if (n_iter > 1) {
          if (ncross != Config.nodes) {
	    /* Incorrect number of crossings: adjust energy  */
	    if ( kkk == 1) {
//	      fprintf(stdout, "Bisection         Energy       Nodes  Discontinuity\n");
	    }
//	    fprintf(stdout, "%5d%25.15e%5d\n", kkk, e, ncross);
	    
	    if (ncross > Config.nodes) {
	      /* Too many crossings: current energy is too high
		 lower the upper bound */
	      eup = e;
	    } else {
	      /* Too few or correct number of crossings: current energy is too low 
		 raise the lower bound */
	      elw = e;
	    }
	    /* New trial value: */
	    e = 0.5 * (eup + elw);
          }
	} else {
 //         fprintf(stdout, "%25.15e%5d%5d\n", e, ncross,Config.nodes);
	}
	
	if ( n_iter == 1 ||  ncross == Config.nodes ) {
	  /*
	    Number of crossings is correct, or energy is fixed:
	    proceed to inward integration 
	    
	    Determination of the wave-function in the last two points 
	    assuming y(Config.mesh+1) = 0
	  */
	  y[Config.mesh] = dx;
	  y[Config.mesh - 1] = (12. - 10.*f[Config.mesh]) * y[Config.mesh] / f[Config.mesh-1];
	  
	  /*	Inward integration */
	  for (i = Config.mesh - 1; i >= icl+1; --i) {
	    y[i-1] = ((12. - 10.*f[i]) * y[i] - f[i+1] * y[i+1]) / f[i-1];
	  }
	  
	  /*	Rescale function to match at the classical turning point (icl) */
	  
	  yicl /= y[icl];
	  for (i = icl; i <= Config.mesh; ++i) {
	    y[i] *= yicl;
	  }
	  
	  /*      normalize on the [-xmax,xmax] segment  */
	  
	  norm = 0.;
	  for (i = 1; i <= Config.mesh; ++i) {
	    norm += y[i]*y[i];
	  }
	  norm = dx * (2.* norm + y[0]*y[0]);
	  norm = sqrt(norm);
	  for (i = 0; i <= Config.mesh; ++i) {
	    y[i] /= norm;
	  }
	  
	  /* 	calculate the discontinuity in the first derivative 
		y'(i;RIGHT) - y'(i;LEFT)         */
	  
	  if (n_iter > 1) {
	    i = icl;
	    djump = (y[i+1] + y[i-1] - (14. - 12.*f[i]) * y[i]) / dx;
//	    fprintf(stdout, "%5d%25.15e%5d%14.8f\n", kkk, e, Config.nodes, djump);
	    if (djump*y[i] > 0.) {
	      /*               Energy is too high --> choose lower energy range */
	      eup = e;
	    } else {
	      /*               Energy is too low --> choose upper energy range */
	      elw = e;
	    }
	    e = 0.5 * (eup + elw);
	  }
	} /* end if (ncross==nodes) */
      } /* end do */
	   
	   /* ---- convergence has been achieved (or it wasn't required) ---- */
      /*
	Calculation of the classical probability density for energy e:
      */
      xmcl = sqrt(2. * e);
      norm = 0.;
      for (i = icl; i <= Config.mesh; ++i) {
	p[i] = 0.;
      }
      for (i = 0; i <= icl - 1; ++i) {
	arg = xmcl*xmcl - x[i]*x[i];
        if ( arg > 0.) 
	  p[i] = 1. / sqrt(arg) / pi;
        else
	  p[i] = 0.;
	norm += dx * 2. * p[i];
      }
      /* The point at x=0 must be counted once: */
      norm -= dx * p[0];
      /* Normalize p(x) so that  Int p(x)dx = 1: */
      for (i = 0; i <= icl - 1; ++i) {
	p[i] /= norm;
      }
      /* lines starting with # ignored by gnuplot */
//      fprintf (out,"#   x       y(x)            y(x)^2       classical p(x)      V\n");
      /* x<0 region: */
      if (hnodes << 1 == Config.nodes)
        parity = +1;
      else
        parity = -1;
	  /*
      for (i = Config.mesh; i >= 1; --i) {
        fprintf(out, "%7.3f%16.8e%16.8e%16.8e%12.6f\n",
		-x[i], parity*y[i], y[i]*y[i], p[i], vpot[i]);
      }*/
      /* x>0 region: */
	  /*
	  for (i = 0; i <= Config.mesh; ++i) {
	      fprintf(out, "%7.3f%16.8e%16.8e%16.8e%12.6f\n",
	      x[i], y[i], y[i]*y[i], p[i], vpot[i]);
      }*/
      /* two blank lines separating blocks of data, useful for gnuplot plotting */
//      fprintf (out,"\n\n");
	
	/* Writing my own output files */
    /* x<0 region: */
   	for (i = Config.mesh; i >= 1; --i) {
       	fprintf(OutFiles.y, "%18.8e ", parity*y[i]);
   	}
   	/* x>0 region: */
   	for (i = 0; i <= Config.mesh; ++i) {
       	fprintf(OutFiles.y, "%18.8e ", y[i]);
   	}
	fprintf(OutFiles.y, "\n");

    	/* writing energy */
	fprintf(OutFiles.e, "%3d%18.8e\n",Config.nodes,e);

	//ADDED
	  calculatedEigenVecNum++;
	}

	closeAllOutFile(OutFiles);
    free(vpot); free(f); free(p); free(y); free(x);
//    fclose(out);
    return 0; // exit(0);
}
