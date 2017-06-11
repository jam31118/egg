#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"

/* Units in this program
 * \tau = \omega * t 
 * (the \tau represent unitless variable and the variable
 * denoted 't' in this progam represent this \tau.
 * The 't' that appears in the above equation represent 
 * the real time in real time unit)
 *
 * \xi = \sqrt{m*\omega / \hbar} * x
 * (the \xi represent unitless that corresponds to 
 * positon variable x with real unit.
 * the 'x' in this script represent \xi and 'x' in
 * the above equation represent the real position variable.)
 *
 * e = E / {\hbar\omega}
 * */

double integralTrapzoid(double *f, double *x, int size) {
	double sum = 0, *pmax = x+size-1;
	for (; x < pmax; x++,f++) {	sum += (*(x+1) - *x) * (*(f+1) + *f); }
	return sum*0.5;
}

double innerprod(double *f, double *g, double *x, int size) {
	double sum = 0, *pmax = x+size-1;
	for (; x < pmax; f++,g++,x++) { 
		sum += (*(x+1) - *x) * ((*(f+1))*(*(g+1)) + (*f)*(*g));
	}
	return sum * 0.5;
}

int main() {

	/* Import CONFIG file */
	Config_t Config;
	char configfile[256] = "CONFIG";
	if(processConfigFile(&Config, configfile)) {
		fprintf(stderr,"[ERROR] During processing CONFIG file\n");
		return 1;
	};
	
	/* Assign config variable to needed variable */
	size_t n; // Number of eigenstates involved
	if (Config.nodesStart == 1) {
		n = Config.nodesEnd;
	} else {
		fprintf(stderr,"[ERROR] Currently unsupported functionality! Please set first index of node to 1\n");
		return 1;
	}
	const size_t N = 2*Config.mesh + 1; // Number of spatial grid points
	const size_t M = Config.tempmesh + 1; // Number of temporal grid points, including \tau = 0
	const size_t cpp_sz = n + n*M + n*N; // Coeff_Phase_Psi vector size
	const double dtau = Config.tmax / Config.tempmesh; // tempmesh is right, not M

	/* Array Declaration and Allocation */
	double *energy = (double *) calloc( n, sizeof(double));
	double *re_Psi = (double *) calloc( N*M, sizeof(double) );
//	double *im_Psi = (double *) calloc( N*M, sizeof(double) );
	double *coeff_cos_re_psi = (double *) calloc( cpp_sz, sizeof(double));
//	double *coeff_sin_im_psi = (double *) calloc( cpp_sz, sizeof(double));
	double *re_psi_trans = (double *) calloc( N*n, sizeof(double)); // form: [.N.][.N.]...[.N.]
//	double *im_psi = (double *) calloc( N*n, sizeof(double)); // form: [.N.][.N.]...[.N.]

	double *psi0 = (double *) calloc( N, sizeof(double));
	double *psi0x = (double *) calloc( N, sizeof(double));

	double *xvec = (double *) calloc( N, sizeof(double));
	double *tau = (double *) calloc( M, sizeof(double));

	/* Array poiners Declaration and Assignment */
	double *re_Psi_p = re_Psi;
	double *re_Psi_p_max = re_Psi + N*M;

	double *coeff_p_head = coeff_cos_re_psi;
	double *cos_p_head = coeff_cos_re_psi + n;
	double *re_psi_p_head = coeff_cos_re_psi + n + n*M;
	
	double *coeff_p_tail = cos_p_head;
	double *cos_p_tail = re_psi_p_head;
	double *re_psi_p_tail = coeff_cos_re_psi + n + n*M + n*N;

	double *coeff_p, *cos_p, *re_psi_p; // *sin_p, *im_psi_p;

	/* Import Eigenfunctions and store into part of coeff_cos_re_psi etc. */
	FILE *feigen = fopen("y.dat","r");
	int idx;
	double *xvec_p, *xvec_p_max = xvec + N;
	double *re_psi_trans_p = re_psi_trans;
	for (xvec_p = xvec; (xvec_p < xvec_p_max) && fscanf(feigen,"%lf",xvec_p); xvec_p++);
	for (idx=0; idx<n; idx++) {
		for (re_psi_p = re_psi_p_head + idx; (re_psi_p < re_psi_p_tail) && fscanf(feigen,"%lf ",re_psi_p); re_psi_p+=n, re_psi_trans_p++) {
			/* saving eigenfunctions in two kinds of data: 
			 * re_psi in coeff_cos_re_psi and re_psi_trans*/
			*re_psi_trans_p = *re_psi_p; 
		}
	}
	/*	
	for (re_psi_p = re_psi_p_head; re_psi_p < re_psi_p_tail; re_psi_p++) {
		fprintf(stderr,"%16.8e",*re_psi_p);
	} fprintf(stderr,"\n");
	*/
	fclose(feigen);

	/* Generate Initial condition (will be relaced by outer module and CONFIG scheme) */
	FILE *f_psi0 = fopen("psi0.dat","w");
//	int indice[3] = {1,2,3}; // ground and 1st exicited
//	double temp = 0;
	for (xvec_p = xvec; (xvec_p < xvec_p_max) && fprintf(f_psi0,"%16.8e",*xvec_p); xvec_p++);
	fprintf(f_psi0,"\n");
	for (idx=0; idx < N; idx++) {
		if ((N/5 < idx) && (idx < N/3)) fprintf(f_psi0,"%16.8e",12.0/N);
	} fprintf(f_psi0,"\n");
	/*
	for (re_psi_trans_p = re_psi_trans; re_psi_trans_p < re_psi_trans + N; re_psi_trans_p++) {
			for (int idx = 0; idx < 3; idx++) {
					temp = (10-indice[idx])/sqrt(14) * (*(re_psi_trans_p+(indice[idx]-1)*N));
			}
			fprintf(f_psi0,"%16.8e",temp);
	}fprintf(f_psi0,"\n");
	*/
	fclose(f_psi0);

	/* Import Initial condition */
	FILE *fpsi0 = fopen("psi0.dat","r");
	double *psi0_p = psi0, *psi0_p_max = psi0 + N;
	double *psi0x_p = psi0x, *psi0x_max = psi0x + N;
	for (; (psi0x_p < psi0x_max) && fscanf(fpsi0,"%lf ",psi0x_p); psi0x_p++);
	if (psi0x_p != psi0x_max) { fprintf(stderr,"[ERROR] Array length mismatch during import\n"); }
	else { fprintf(stderr, "[ LOG ] Array length was matched during import\n"); }
	for (; (psi0_p < psi0_p_max) && fscanf(fpsi0,"%lf ",psi0_p); psi0_p++ );
	if (psi0_p != psi0_p_max) { fprintf(stderr,"[ERROR] Array length mismatch during import\n"); }
	else { fprintf(stderr, "[ LOG ] Array length was matched during import\n"); }
	/*
	for (psi0x_p = psi0x; psi0x_p < psi0x_max; psi0x_p++) {
		fprintf(stderr,"%16.8e ",*psi0x_p);
	} fprintf(stderr,"\n");
	for (psi0_p = psi0; psi0_p < psi0_p_max; psi0_p++) {
		fprintf(stderr,"%16.8e ",*psi0_p);
	} fprintf(stderr,"\n");
	*/
	fclose(fpsi0);
	

	/* Import eigen energy data, e_i, i \in {1,...,n} */
	FILE *fp_energy = fopen("energy.dat","r");
	double *energy_p, *energy_p_max = energy + n;
	for (energy_p = energy; (energy_p < energy_p_max) && fscanf(fp_energy,"%*d %lf",energy_p); energy_p++);
	if (energy_p != energy_p_max) {fprintf(stderr,"[ERROR] Array length mismatch during import\n");}
	else {fprintf(stderr,"[ LOG ] Array length was matched during energy value import\n");}
	/*
	for (energy_p =energy; energy_p < energy_p_max; energy_p++) {fprintf(stderr,"%16.8e",*energy_p);}
	*/
	fclose(fp_energy);

	/* Assigning \tau variable with equal timestep */
	double *tau_p, *tau_p_max = tau + M, tau_val = 0;
	for(tau_p=tau; tau_p < tau_p_max; tau_p++, tau_val += dtau) {*tau_p = tau_val;}

	/* Calculate cos(e_i * tau_j) and put it to cos part in coeff_cos_re_psi */
	cos_p = cos_p_head;
	for (tau_p = tau; tau_p < tau_p_max; tau_p++) {
		for (energy_p = energy; energy_p < energy_p_max; energy_p++, cos_p++) {
			*cos_p = cos(*energy_p * (*tau_p)); 
			// [THINK] *tau_p may invoke LW instr. is it good?
			// Can it be optimized further?
			//fprintf(stderr,"%16.8e",*cos_p);
		}
	}
	//for(cos_p=cos_p_head; cos_p < cos_p_tail; cos_p++) {fprintf(stderr,"%16.8e",*cos_p);}
	/* (170520 sin, or imaginary, version should also be implemented )*/

	/* Calculate coefficient and put them into coeff part in coeff_cos_re_psi */
	re_psi_trans_p = re_psi_trans;
	//double *tmpN = (double *) calloc(N,sizeof(double));
	//double *tmpN_p = tmpN, tmpN_p_max = tmpN + N;
	for(coeff_p = coeff_p_head; coeff_p < coeff_p_tail; coeff_p++, re_psi_trans_p+=N) {
		//for(tmpN_p = tmpN; tmpN_p < tmpN_p_max; tmpN++) {}
		//*coeff_p = integralTrapzoid(re_psi_trans_p .. ,,N);
		*coeff_p = innerprod(psi0, re_psi_trans_p, xvec, N); 
	}
	/* (170520 sin, or imaginary, coefficient should also be implemented )*/
		
	/* Martrix multiplication for propagation */
	for (cos_p = cos_p_head; cos_p < cos_p_tail; cos_p += n) {
		for (re_psi_p = re_psi_p_head; re_psi_p < re_psi_p_tail; )
		// there should be no re_psi_p increment since it is in the lower loop.
		{
			for (coeff_p = coeff_p_head; coeff_p < coeff_p_tail; coeff_p++) {
				*re_Psi_p += (*coeff_p) * (*cos_p) * (*re_psi_p);
				cos_p++; // ADDI rather than ADDing large index
				re_psi_p++; // ADDI rather than ADDing large index
			//	fprintf(stdout,"[ LOG ] Running with coeff_p == %p\n",coeff_p);
			}
			cos_p -= n;
			re_Psi_p++;
		}
		//fprintf(stderr,"[ LOG ] Left data to be processed: %li\n",cos_p_tail - cos_p);
	}


	/* Print the result re_Psi to file */
	FILE *f_tdse = fopen(Config.filename.tdseout,"w");
	for (xvec_p = xvec; xvec_p < xvec_p_max; xvec_p++) {
		fprintf(f_tdse,"%18.8e",*xvec_p);
	} fprintf(f_tdse,"\n");
	double *tmp = re_Psi;
	for (re_Psi_p = re_Psi; re_Psi_p < re_Psi_p_max; re_Psi_p += N) {
		for (tmp = re_Psi_p; tmp < re_Psi_p + N; tmp++) {
			fprintf(f_tdse,"%18.8e",*tmp);
		} fprintf(f_tdse,"\n");
	}
	fclose(f_tdse);
	
}
