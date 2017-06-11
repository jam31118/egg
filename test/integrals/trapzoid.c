#include <stdio.h>

void printArr(double *x, int size) {
	for (int i=0; i<size; i++) {
		printf("x == %.2f\n",*(x+i));
	}
}

double integralTrapzoid1D(double *x, double *f, int size) {
	double *x_p, sum = 0;
	double *f_p = f;
	printf("sum in func == %.2f\n",sum);
	printArr(f,size);
	for (x_p = x; x_p < x+size-1; x_p++) {
		sum += (*(x_p+1) - *x_p) * (*(f_p+1) + *f_p);
		printf("f == %.2f\tf+1 == %.2f\n",*f_p,*(f_p+1));
		printf("sum in func == %.2f\n",sum);
		f_p++;
	}
	return sum*0.5;
}

double integralTrapzoid(double *x, double *f, int size) {
	double sum = 0, *pmax = x+size-1;
	for (; x < pmax; x++,f++) {	sum += (*(x+1) - *x) * (*(f+1) + *f); }
	return sum*0.5;
}
void chptr(int *a) {
	printf("a == %p\n",a);
	a++;
	printf("a == %p\n",a);
}

int main() {
	double x[3] = {0,1,2};
	double f[3] = {1,1,2};

	printArr(x,3);
	printArr(f,3);

	//double sum = integralTrapzoid1D(x,f,3);
	double sum = integralTrapzoid(x,f,3);
	printf("sum == %.2f\n",sum);
	/*
	int xx = 3;
	int *xxp = &xx;
	printf("xxp == %p\n",xxp);
	chptr(xxp);
	printf("xxp == %p\n",xxp);
	*/
	
}
