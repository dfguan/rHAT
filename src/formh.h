
#include <stdint.h>

#define PATH_LEN 1024

typedef struct options {
	uint32_t 	len_sed;
	char 		refpath[PATH_LEN];
	char 		hashdir[PATH_LEN];

}opts;


class Form {
	
	opts *opt;
public:
	Form(opts *opt);
	int usage();
	int opt_parse(int argc, char *argv[], opts* opt);
};