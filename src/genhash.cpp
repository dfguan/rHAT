#include <stdio.h>
#include "formh.h"
#include "whash.h"
#include "readfl.h"
#include <iostream>
#include <time.h>
#include <cstdlib>
using namespace std;

const std::string getCurrentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "[INFO] %Y-%m-%dT%X", &tstruct);

    return buf;
}

int main(int argc, char *argv[])
{
	
	opts *opt = new opts;
	Form fm(opt);
	if (fm.opt_parse(argc,argv,opt)!=1)
		exit(1);

	fprintf(stderr,"%s rHAT-indexer started\n",getCurrentDateTime().c_str());
	
	uint32_t len_genome = 0;
	
	read_file rdfl;
	char 	*genome = rdfl.read_ref(opt->refpath,&len_genome);

	Hash hashh;
	hashh.write_hashfile(opt->hashdir,genome,len_genome,opt->len_sed);

	
	if ( NULL != opt ) delete opt;
	fprintf(stderr,"%s rHAT-indexer ended\n",getCurrentDateTime().c_str());
	return 0;

}
