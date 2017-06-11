#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "./config.h"

int processConfigFile(Config_t *config_p, const char *filename) {
	fprintf(stderr, "[ LOG ] Starting config processing\n");

	FILE *fp = fopen(filename,"r");
	if (fp == NULL) {
		fprintf(stderr,ERR_FILEOPEN);
		return 1;
	}
	/* Read CONFIG file */
	char line[bufsize];
	char *key = NULL, *value = NULL;
	while ( fgets(line,bufsize,fp) != NULL ) {
		//fprintf(stderr,"[ LOG ] reading CONFIG file line by line\n");
		/* Passes blank line */
		if (!strcmp(line,"\n")) { continue; }

		/* Parsing KEY and VALUE */
		key = strtok(line,delim);
		if (!key) {
			fprintf(stderr,ERR_TOKENIZE);
			fprintf(stderr,"[ERROR] Cannot get KEY\n");
			return 1;
		}
		value = strtok(NULL,delim);
		if (!value) {
			fprintf(stderr,ERR_TOKENIZE);
			fprintf(stderr,"[ERROR] Cannot get VALUE\n");
			return 1;
		}
		/* Tries to store the parsed value into Config struct. 
		 * If there's no matching key among predefined config
		 * or other problem arised during storing,
		 * stop this program with error log. */
		if (storeConfig(config_p,key,value)) {
			fprintf(stderr,"[ERROR] Abnormal config storing\n");
			return 1;
		}
	}
	fclose(fp);
	
	fprintf(stderr,"[ LOG ] CONFIG file processed completed\n");
	return 0;
}

int storeConfig(Config_t *config_p, char *key, char *value) {
	char *endptr = "";
	if (!strcmp(key,"xmax")) {
		config_p->xmax = strtod(value,&endptr);
	} else if (!strcmp(key,"tmax")) {
		config_p->tmax = strtod(value,&endptr);
	} else if (!strcmp(key,"mesh")) {
		config_p->mesh = strtol(value,&endptr,10);
	} else if (!strcmp(key,"tempmesh")) {
		config_p->tempmesh = strtol(value,&endptr,10);	
	} else if (!strcmp(key,"nodes")) {
		/* 170519 implementation 
		 * the value is assumed to be the format:
		 * 'n:m' where n, m are indice of nodes(<->eigenstate)
		 * and it calculate m-n+1 eigenstate with #ofNodes
		 * of n, n+1, ..., m */
		/* check if it is in right format 'n:m' */
		char *tok = strtok(value,":");
		if (tok == NULL) {
			fprintf(stderr,"[ERROR] Out of format!\n");
			return 1;
		}
		/* get first index of 'n:m', namely, n */
		config_p->nodesStart = strtol(tok,&endptr,10);
		/* get last index of 'n:m', namely, m */
		tok = strtok(NULL,":");
		config_p->nodesEnd = strtol(tok,&endptr,10);

        // 170514 should be replaced
		// 170519 replaced
		// config_p->nodes = strtol(value,&endptr,10);
		
        /* 170514 new implementation: 
         * int* in config will be assigned a pointer to
         * integer array of nodesIndices */
        /* ... */
	} else if (!strcmp(key,"fileout")) {
		strcpy(config_p->fileout, value);
	} else if (!strcmp(key,"tdseout")) {
		strcpy(config_p->filename.tdseout, value);
	} else if (!strcmp(key,"e")) {
		config_p->e = strtod(value,&endptr);
	} else if (!strcmp(key,"videoname")) {
		fprintf(stderr,"[ LOG ] videoname has been enterd in CONFIG. This would be implemented later\n");
	}
	else {
		fprintf(stderr,"[ERROR] no matching key\n");
		return 1;
	}

	/* Check if the value is valid integer string 
	 * If *endpth == '\0' then the value is valid */
	if (*endptr != '\0') {
		fprintf(stderr,"[ERROR] something wrong in parsing\n");
		return 1;
	}

	return 0;
}

