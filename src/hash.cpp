#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "hash.h"

uint8_t trans[128]= {
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4
	};

int comparator(const void *p, const void *q)
{
	seed 	*t = (seed *)p;
	seed 	*f = (seed *)q;
	if (t->sed_value > f->sed_value)
		return 1;
	else {
		if(t->sed_value == f->sed_value) {
			if (t->sed_bkt_pos > f->sed_bkt_pos) {
				return 1;
			} else {
				if(t->sed_bkt_pos == f->sed_bkt_pos)
					return 0;
				else
					return -1;
			}

		}
		else
			return -1;
	}
}

int transfer(char *genome,uint32_t *l2r,uint32_t len_sed)
{

	*l2r = 0;
	uint32_t temp;
	for(uint32_t i = 0;i<len_sed;i++) {
		temp = trans[genome[i]];
		*l2r = *l2r << 2;
		*l2r = *l2r | temp;
	}
	return 0;
}

uint32_t statsed(seed *repsed,uint32_t len)
{	
	uint32_t count = 1;
	uint32_t temp_sed_value = repsed[0].sed_value;
	uint32_t temp_sed_bkt_pos = repsed[0].sed_bkt_pos;
	for (uint32_t i=1;i<len;++i) {
		if (repsed[i].sed_value == temp_sed_value && repsed[i].sed_bkt_pos == temp_sed_bkt_pos )
			continue;
		else {
			++count;
			temp_sed_value = repsed[i].sed_value;
			temp_sed_bkt_pos = repsed[i].sed_bkt_pos;			
		}
	}
	return count;

}


Hashtab *Hash::load_hashfile(char *path,uint32_t len_genome, uint32_t len_sed)
{
	//char num[4] = {0};
	uint32_t amount = 1 << (len_sed<<1);

	Hashtab * hashtab = new Hashtab;//may be something wrong;
	if( NULL == hashtab) {
		fprintf(stderr, "Failed when applying new space for hash table");
		exit(1);
	}

	hashtab->pointer = new uint32_t[amount + 1];//may be something wrong;

	if ( NULL == hashtab->pointer ) {
		fprintf(stderr, "Failed when applying new space for hash table");
		exit(1);

	}

	for(uint32_t i=0;i<=amount;++i) {
		hashtab->pointer[i] = 0;
	}
	//error ask for too much memory

	/*hashtab->seq_bkt = new char [3][len_genome - len_sed + 1];//may be something wrong;
	for (uint32_t i=0; i<=len_genome - len_sed; ++i) {
		hashtab->seq_bkt[i] = new char[3];
	}
	*/
	hashtab->seq_bkt = new uint32_t [len_genome-len_sed + 1];

	if ( NULL == hashtab->seq_bkt ) {
		fprintf(stderr, "Failed when applying new space for hash table");
		exit(1);
	}

	char filepath[256] = {0};

	FILE *fp_hash;

	strcpy(filepath,path);
	strcat(filepath,"/");
	strcat(filepath,"hash");
	

	if ((fp_hash = fopen(filepath,"rb")) == NULL) {
		fprintf(stderr,"Failed to open hash file, now exit");
		exit(1);
	}
	//something may be wrong here cause sizeof(uint32_t)*amount may be out of range	
	fread(hashtab->pointer,sizeof(uint32_t),amount+1,fp_hash);
	fread(hashtab->seq_bkt,sizeof(uint32_t),len_genome-len_sed+1,fp_hash);
	fclose(fp_hash);
	return hashtab;
}
