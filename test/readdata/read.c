#include <stdio.h>
#include <stdlib.h>


int main() {

	/* Read from file .. N etc.*/
	int N = 50; 
	
	/* Declare an array */
	double *aa = (double *) calloc(2*N+1,sizeof(double));
	double *xx = (double *) calloc(2*N+1,sizeof(double));
	double *aamax = aa + 2*N+1;
	double *aa_p = aa;
	double *xx_p = xx;

	/* Import data file */
	FILE *fp = fopen("./y.dat","r");

	/*
	aa_p = aa;
	while (fscanf(fp,"%lf ",aa_p)) {

		aa_p++;
	}*/
	/*
	for (aa_p = aa; (aa_p < aamax) && fscanf(fp,"%lf ",aa_p); aa_p++) {
		//fscanf(fp,"%lf ",aa_p);
	}
	if (aa_p != aamax) { fprintf(stderr, "[ERROR] Mismatching array length in data import"); }
	*/
	int num =0; 
	char c = 0; // let's check whether it is '\n' or not
	for (aa_p = aa; (aa_p <aamax) && fscanf(fp,"%lf ",aa_p); aa_p++) {
	//	fprintf(stderr,"%p\t%lf\tnum=%d\n",aa_p,*aa_p,num);
	//	if (aa_p > aamax+10) break;
		//fscanf(fp,"%lf ",aa_p);
	}

	if (num = fscanf(fp,"%c",&c)) {
		if (c=='\n') { printf("oh!,");}
		else {printf("no \\n, instead it is == %c\twith num == %d\n",c,num);}
	} else {printf("nothing?\n"); }
	if (num = fscanf(fp,"%c",&c)) {
		if (c=='\n') { printf("oh!,");}
		else {printf("no \\n, instead it is == %c\twith num == %d\n",c,num);}
	} else {printf("nothing?\n"); }
	if (num = fscanf(fp,"%c",&c)) {
		if (c=='\n') { printf("oh!,");}
		else {printf("no \\n, instead it is == %c\twith num == %d\n",c,num);}
	} else {printf("nothing?\n"); }

//	if (aa_p != aamax) { fprintf(stderr, "[ERROR] Mismatching array length in data import"); }
	fclose(fp);

	/* Print the result */ 
	for (aa_p = aa; aa_p < aamax; aa_p++) {
		fprintf(stdout,"%16.8e ",*aa_p);
	}
	fprintf(stdout,"\n");
}
