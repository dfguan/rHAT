
#ifndef ALIGNER_H_
#define ALIGNER_H_
#define LEN_BASES 1024

#define LEN 100


#include "form.h"
#include "hash.h"

#include "graph.h"

#include "kseq.h"

#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

#include <iostream>
#include <cstdio>
#include <cstring>
#include <queue>
using namespace std;

typedef struct bucket
{
	uint32_t 		array_seq;// actually how many bucket should there be;
	uint32_t 		sequence_num_bkt;// also don't know how many
	bool operator<(const bucket& r) const
	{
		return sequence_num_bkt > r.sequence_num_bkt;
	}

}bkt;

typedef struct bucket2
{
	uint16_t	hit_times;
	uint32_t 	seq_num;
	bool		isrc;
	bool operator<(const bucket2 & r) const
	{
		return (hit_times < r.hit_times) ||(hit_times == r.hit_times && seq_num < r.seq_num);
	}
	bucket2 operator=(const bucket2 & r) 
	{	
		hit_times = r.hit_times;
		seq_num = r.seq_num;
		isrc = r.isrc;
		return *this;
	}
}bkt2;




class Aligner {
	
	bkt2 		preserved[20];//should be a defined number write
	uint32_t 	pos[20];//write	
	uint8_t 	chrIndex[20];//write
	int 		countChr;//readable  
	char 		*genome;//readable
	uint32_t 	len_genome;//readable
	Hashtab 	*hashtab;//readable
	char 		**ChrName;//readable
	uint32_t 	Start_pos[100];//readable
	opts 		*opt;//readable
	int8_t 		mat[25];
public:
	Aligner(opts *opt);
	~Aligner();
	void Runtask();
private:
	uint32_t 		*sed_rec;//len_bases -14] ;//use malloc();
	uint16_t 		*sed_hit_times;//use malloc
	uint16_t 		*unused_bkt;
	char 			*RCRead;
	Sam_Rec 		*sams;	
	Sam_Rec 		*svsams;
	int 			applyNonSV(kseq_t *trunk,  RHashtable *rhashtab, RHashtable *rrhashtab);
	void 			applySV(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab);
	int  			conductAlign(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab, 
					std::priority_queue <bkt2> &cansHeap) ;
	int 			conductAlign(kseq_t *trunk, char *read, char *rcRead, int lenRead, RHashtable *rhashtab, RHashtable *rrhashtab
					,std::priority_queue <bkt2> &cansHeap, Sam_Rec *svsamsp);
	void 			proCans(char *read, uint32_t lenRead, bool isRC, std::priority_queue<bkt2> &cansHeap);

	uint32_t 		gen_sed(char *bases, uint32_t lenRead);
	
	void 			sort_sed(uint32_t usedseed);
	void  			stat_sed(uint32_t usedseed);
	int 			RevComRead(char *read,int len_read);
	bool 			IsQualifiedRead(char *read,int len);

	//sv

	int 			Qualified(Sam_Rec *key,Sam_Rec *set, int start, int end, bool *issubstr);
	int 			connect(Sam_Rec *rec1, Sam_Rec *rec2, kseq_t *trunk);
	int 			produceSAM(int countbulks,int *sam4bulk, kseq_t *trunk, uint32_t *len);
};
extern 	uint8_t rev[];
#endif
