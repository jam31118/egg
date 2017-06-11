#include <stdio.h>
#include <string.h>

int main() {
    char str[40] = "1,2,3,5:6,9:23";
    printf("%s\n",str);
    char *tok = NULL;
    tok = strtok(str,",");
    printf("tok: %s\n",tok);
    if (!tok) {
        fprintf(stderr,"[ERROR] cannot tokenize\n");
        return 1;
    }
//    char *tok1 = NULL;
    while (tok = strtok(NULL,",")) {
        printf("tok: %s\n",tok);
        if (strchr(tok,':')) {
            fprintf(stderr,"[ LOG ] Found ':'\n");
        }
    }
    return 0;
}
