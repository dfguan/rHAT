
#include "form.h"
#include "desc.h"

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

#define RH_SEED_LEN 8
#define WAITING_LEN 203



char *const short_options = "t:n:hl:m:g:o:e:c:d:";
struct option long_options[] = {
    { "threads",     1,   NULL,    't'   },
    { "num",     1,   NULL,    'n' },
    { "help",     0,   NULL,    'h'   },
    { "seed_length",     1,   NULL,    'l'   },
    { "hit_max",		1,NULL,	'm'},
    //{ "auto_load", 0, NULL, 'a'},
    {"longest", 1,  NULL, 'g'},
    {"gapopen", 1,  NULL,'o'},
    {"gapextended", 1,  NULL,'e'},
    {"match",   1,  NULL,'c'},
    {"mismatch",    1,  NULL,'d'},
    { 0,     0,   0,    0   }
};

Form::Form(opts *opt)
{

    opt->len_sed = 13;
    opt->canN = 5;
    opt->hit_limit = -1;
    //opt->usecigar = true;
    opt->rh_seed_len = RH_SEED_LEN;
    opt->waitingLen = WAITING_LEN;
    opt->thread = 1000;
    opt->len_limit = 50000;
    opt->gapopen = 2;
    opt->gapextend = 1;
    opt->mismatch = 5;
    opt->match = 5;

    //opt->autoload = false;

}

int Form::usage()
{

    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   rHAT-mapper\n"); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
    fprintf(stderr, "Usage:     rHAT-mapper [Options] <IndexDir> <Read> <Reference>\n\n"); 
    fprintf(stderr, "Options:   -h, --help                   help\n"); 
    fprintf(stderr, "           -t, --threads       <int>    thread\n"); 
    fprintf(stderr, "           -n, --num           <int>    candidate number [5]\n"); 
    fprintf(stderr, "           -m, --hit_max       <int>    max hit times of a seed [65,535]\n"); 
    fprintf(stderr, "           -l, --seed_length   <int>    seed length of hash index [13]\n"); 
    fprintf(stderr, "           -g, --longest       <int>    longest length of read [50,000]\n");
    fprintf(stderr, "           -o, --gapopen       <int>    gapopen penalty of ksw [2]\n");
    fprintf(stderr, "           -e, --gapextended   <int>    gapextended penalty of ksw [1]\n");
    fprintf(stderr, "           -c, --match         <int>    match score of ksw [1]\n");
    fprintf(stderr, "           -d, --mismatch      <int>    mismatch score of ksw [5]\n");
    //fprintf(stderr, "           -a, --auto_load              load hash table from hash file without produce hash file");
    //fprintf(stderr, "           -c, --write_cigar            print cigar in XA fields [False]\n"); 
    fprintf(stderr, "\n"); 
    return 0;
}

int Form::opt_parse(int argc, char *argv[], opts* opt)
{
	int c; 
    int option_index=0;
    if (argc == 1) return usage();
    while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0){
        switch(c){
            case 't':
                opt->thread = atoi(optarg);
                break;
            case 'n':
                opt->canN = atoi(optarg);
                break;
            case 'h':
                return usage();
                break;
            case 'm':
                opt->hit_limit = atoi(optarg);
                break;
            case 'l':
                opt->len_sed = atoi(optarg);
                break;
            case 'g':
                opt->len_limit = atoi(optarg);
                break;
            case 'o':
                opt->gapopen = atoi(optarg);
                break;
            case 'e':
                opt->gapextend = atoi(optarg);
                break;
            case 'c':
                opt->match = atoi(optarg);
                break;
            case 'd':
                opt->mismatch = atoi(optarg);
                break;
            default:
                fprintf(stderr,"not proper parameters\n");
                return usage();
              
        }
    
    }
    if(optind + 3 != argc){
        fprintf(stderr, "[opt_parse]: index directory, read file and reference file can't be omited!\n"); 
        return 0; 
    }
    opt->argv = argv;
    opt->argc = argc;
 	strncpy(opt->hashdir, argv[optind++],sizeof(opt->hashdir));
 	strncpy(opt->readpath,argv[optind++],sizeof(opt->readpath));
 	strncpy(opt->refpath,argv[optind++],sizeof(opt->refpath));

    return 1;  

}
