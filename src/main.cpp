/*
	Description: first parameter: lenght of seed 
	second parameter: choice of producing hash table
	third paramter: number of reference area
	forth paramter: whether writing cigar
	fifth: the path of comparision file;
*/
#include "aligner.h"

#include <time.h>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
using namespace std;

#define PRINT_LOG 

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

#ifdef PRINT_LOG
	fprintf(stderr,"%s [rHAT-mapper] started\n",getCurrentDateTime().c_str());
#endif
	
	opts *opt = new opts;
	Form fm(opt);
	if (fm.opt_parse(argc,argv,opt)!=1)
		exit(1);
	
	Aligner alig(opt);
	alig.Runtask();	

	if ( NULL != opt ) delete opt;

#ifdef PRINT_LOG
	fprintf(stderr,"%s [rHAT-mapper] ended\n",getCurrentDateTime().c_str());
	//fclose(stderr);
#endif	
	return 0;

}

