#include <stdio.h>

int main() {

	int i,j;
	j = 0;
	int a[3] = {-3,-2,-1};
	for (i=0; i<3; i++) {
		printf("a[%d]=%d\n",i,a[i]);
		j++;
	}
	printf("i:%d\tj:%d\n",i,j);
	printf("a[%d]=%d\n",i,a[i]);

	for (i=0; i<3; ) {
		printf("i=%d\n",i);
		i++;
	}
}
