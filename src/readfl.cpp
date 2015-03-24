#include "readfl.h"
#include <iostream>
#include <fstream>
using namespace std;

#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

char changeRef[128] = {
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','A','N','C','N','N','N','G','N','N','N','N','N','N','N','N',
    'N','N','N','N','T','N','N','N','N','N','N','N','N','N','N','N',
    'N','A','N','C','N','N','N','G','N','N','N','N','N','N','N','N',
    'N','N','N','N','T','N','N','N','N','N','N','N','N','N','N','N'
};

void read_file::preDealRef(char *ref, uint32_t len_ref) 
{
    for (uint32_t i=0;i<len_ref;++i) {
        ref[i] = changeRef[ref[i]];
    }

}

void combineStr(char *str1,char *str2,uint32_t len) 
{
    for (uint32_t i=0;i<len;++i) {
        str1[i] = str2[i];

    }
    str1[len] = '\0';

}

char *read_file::read_ref(char *path,uint32_t *len_genome, uint32_t *Start_pos, char **chrName,int *count_chrl)
{
        char *genome = NULL;
    //get the size of file 
        ifstream fl;
        fl.open(path);
        if (fl.fail()) {
            fprintf(stderr,"Failed to open reference file, now exit\n");
            exit(1);
        } 
        streampos begin,end;
        begin = fl.tellg();
        fl.seekg(0,ios::end);
        end = fl.tellg();
        unsigned flsize = (unsigned)(end - begin);
        fl.close();
    //use kseq to read ref.  
        genome = new char[flsize];
        gzFile fp;
        kseq_t *trunk; 
        fp = gzopen(path, "r");
        trunk = kseq_init(fp);
        //uint32_t size = -1;
        
        int index = 1;
        Start_pos[0] = 0;
        if ( NULL == genome) {
           fprintf(stderr,"Failed to applying for new space, now exit");
            exit(1);
        }
        while(kseq_read(trunk) >= 0) {
            //genome = new char[trunk->seq.l+1];
            //*len_genome += trunk->seq.l;
            combineStr(genome+*len_genome,trunk->seq.s,trunk->seq.l);
            *len_genome += trunk->seq.l;
            Start_pos[index] = *len_genome;
            strcpy(chrName[index++],trunk->name.s);
        }
        *count_chrl = index;
        preDealRef(genome,*len_genome);
        kseq_destroy(trunk);
        gzclose(fp);
        //fprintf(stderr,"%u",*len_genome);

        return genome;
}

char *read_file::read_ref(char *path,uint32_t *len_genome)
{
        char *genome = NULL;
    //get the size of file 
        ifstream fl;
        fl.open(path);
        if (fl.fail()) {
            fprintf(stderr,"Failed to open reference file, now exit\n");
            exit(1);
        } 
        streampos begin,end;
        begin = fl.tellg();
        fl.seekg(0,ios::end);
        end = fl.tellg();
        unsigned flsize = (unsigned)(end - begin);
        fl.close();
    //use kseq to read ref.  
        genome = new char[flsize];
        gzFile fp;
        kseq_t *trunk; 
        fp = gzopen(path, "r");
        trunk = kseq_init(fp);
        //uint32_t size = -1;
        
        //int index = 1;
        //Start_pos[0] = 0;
        if ( NULL == genome) {
             fprintf(stderr,"Failed to applying for new space, now exit");
            exit(1);
        }
        while(kseq_read(trunk) >= 0) {
            //genome = new char[trunk->seq.l+1];
            //*len_genome += trunk->seq.l;
            combineStr(genome+*len_genome,trunk->seq.s,trunk->seq.l);
            *len_genome += trunk->seq.l;
            //Start_pos[index] = *len_genome;
            //strcpy(chrName[index++],trunk->name.s);
        }
        //*count_chrl = index;
        preDealRef(genome,*len_genome);
        kseq_destroy(trunk);
        gzclose(fp);
        //fprintf(stderr,"%u",*len_genome);

        return genome;
}
