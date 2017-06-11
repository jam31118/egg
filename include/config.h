#ifndef _CONFIG_H_
#define _CONFIG_H_

#define FILENAMEMAX 256

/* Define Constants */
static const unsigned int bufsize = 256;
static const char delim[] = "=\n";

/* Define Configuration structure 
 * storeConfig() should be sync. with this struct. */
typedef struct Config_t {
	double xmax, tmax;
	int mesh, tempmesh, nodes;
	int nodesStart, nodesEnd;
	char fileout[FILENAMEMAX];
	struct {
		char tdseout[FILENAMEMAX];
	} filename;
	double e;
} Config_t;

/* Define ERROR messages */
static const char ERR_FILEOPEN[] = "[ERROR] Cannot open file properly\n";
static const char ERR_TOKENIZE[] = "[ERROR] Abnormal tokenization\n";

/* Function Declaration */
int processConfigFile(Config_t *config_p, const char *filename);
int storeConfig(Config_t *config_p, char *key, char *value);

#endif
