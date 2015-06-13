
#include "formh.h"
#include "desc.h"


#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>




char *const short_options = "k:h";
struct option long_options[] = {
    //{ "threads",     1,   NULL,    't'   },
    //{ "num",     1,   NULL,    'n' },
    { "kmer-size",     1,   NULL,    'l'   },
    { "help",     0,   NULL,    'h'   },
    
    //{ "hit_max",		1,NULL,	'm'},
    //{ "auto_load", 0, NULL, 'a'},
    { 0,     0,   0,    0   }
};

Form::Form(opts *opt)
{

    opt->len_sed = 13;

}

int Form::usage()
{

    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   rHAT-indexer\n"); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
    fprintf(stderr, "Usage:     rHAT-indexer [Options] <HashIndexDir> <Reference>\n\n"); 

    fprintf(stderr, "<HashIndexDir>         The directory storing RHT index\n");
    fprintf(stderr, "<Reference>            Sequence of reference genome, in FASTA format\n\n");
    
    fprintf(stderr, "Options:   -k, --kmer-size        <int>           the size of the k-mers extracted from reference genome for indexing [13]\n"); 
    fprintf(stderr, "           -h, --help                             help\n");
    //fprintf(stderr, "           -t, --threads       <int>    thread\n"); 
    //fprintf(stderr, "           -n, --num           <int>    candidate number [5]\n"); 
    //fprintf(stderr, "           -m, --hit_max       <int>    max hit times of a seed [1000]\n"); 
    //fprintf(stderr, "           -l, --seed_length   <int>    seed length of hash index [13]\n"); 
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
            case 'h':
                return usage();
                break;
            case 'k':
                opt->len_sed = atoi(optarg);
                break;
            default:
                fprintf(stderr,"not proper parameters\n");
                return usage();
              
        }
    
    }
    if(optind + 2 != argc){
        fprintf(stderr, "[opt_parse]: index directory and reference file can't be omited!\n"); 
        return 0; 
    }

 	strncpy(opt->hashdir, argv[optind++],sizeof(opt->hashdir));

 	strncpy(opt->refpath,argv[optind++],sizeof(opt->refpath));

    return 1;  

}
