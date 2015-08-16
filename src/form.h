
#include "desc.h"

#include <stdint.h>

#define PATH_LEN 1024

typedef struct options {
	uint32_t 	len_sed;
	uint32_t 	canN;
	char 		readpath[PATH_LEN];
	char 		refpath[PATH_LEN];
	char 		hashdir[PATH_LEN];
	uint16_t 	hit_limit;
	uint32_t 	rh_seed_len;
	uint32_t 	len_limit;
	uint32_t 	waitingLen;
	int 		thread;
	int 		argc;
	char 		**argv;
	//int 		localKmer;
	int 		gapopen;
	int			gapextend;
	int			match;
	int 		mismatch;
	
}opts;


class Form {
	
	opts *opt;
public:
	Form(opts *opt);
	int usage();
	int opt_parse(int argc, char *argv[], opts* opt);
};
