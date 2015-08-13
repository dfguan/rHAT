
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

typedef struct {
	char	*qual;
	int 	qual_len;

	char   	*seq;
	char 	*rseq;
	int 	seq_len;

	char 	*name;
	int 	name_len;
}rhat_seq;

typedef struct {
	uint32_t		*sed_rec;
	uint16_t 		*sed_hit_times;
	uint16_t		*unused_bkt;

	RHashtable 		*rhashtab;
	RHashtable 		*rrhashtab;

}aux_var;

class Aligner {
	bkt2 		preserved[20];//should be a defined number write
	uint32_t 	pos[20];//write	
	uint8_t 	chrIndex[20];//write
	int 		countChr;//readable  
	char 		*genome;//readable
	char 		*genome_e;
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
	
	int 			applyNonSV(kseq_t *trunk,  RHashtable *rhashtab, RHashtable *rrhashtab, Sam_Rec *_sams, 
	uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt);

	int 			applySV(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab, SvSam_Rec *_svsams, 
	uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt);
private:
	uint32_t 		*sed_rec;//len_bases -14] ;//use malloc();
	uint16_t 		*sed_hit_times;//use malloc
	uint16_t 		*unused_bkt;
	
	//uint32_t 		sv_interval;
	

	int 			conductAlign(kseq_t *trunk, char *read, char *rcRead, int lenRead, RHashtable *_rhashtab, RHashtable *_rrhashtab
	,std::priority_queue <bkt2> &cansHeap, SvSam_Rec *_svsamsp); 

	int 			conductAlign(kseq_t *trunk,std::priority_queue <bkt2> &cansHeap, RHashtable *_rhashtable, RHashtable *_rrhashtable, 
		Sam_Rec *_sams); 
	
	void 			proCans(char *read, uint32_t lenRead, bool isRC,  std::priority_queue<bkt2> &cansHeap, 
	uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt) ;

	uint32_t 		gen_sed(char *bases, uint32_t len_bases, uint32_t *_sed_rec);
	
	void 			sort_sed(uint32_t *_sed_rec,uint32_t usedseed);

	void  			stat_sed(uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt, uint32_t usedseed);
	
	int 			revComRead(char *read, char *rcRead, int len_read);
	
	bool 			isQualifiedRead(char *read,int len);

	
	//sv
	int 			rhat_seq_read(kstream_t *_fp, kseq_t *_seqs, int n_needed);
	int 			qualified(SvSam_Rec *key,SvSam_Rec *set, int start, int end);
	int 			connect(SvSam_Rec *rec1, SvSam_Rec *rec2, kseq_t *trunk);
	int 			produceSAM(SvSam_Rec *_svsams , int countbulks,int *sam4bulk, kseq_t *trunk, uint32_t *len);
	int 			mergeBoundaryCigar(string & cigar1, string & cigar2, uint32_t *cigar, int n_cigar, const char *correspondTable);
	//output
	int 			outputSam(kseq_t *_seqs, Sam_Rec *_sams, SvSam_Rec **_svsams, uint16_t *_sam_details, int _n_seqs);
};
typedef struct aux
{
	int 		tid;
	Aligner 	*aln;
	int 		n_seqs;
	aux_var 	com_var;
	kseq_t 		*seqs;
	Sam_Rec		*sams;
	SvSam_Rec	**svsams;
	opts 		*opt;
	uint16_t 	*sam_details;// H:6 samamount L:2 isSAM, bits set by limit of read length might be changed in future
}thread_aux;
thread_aux *thread_initiate(int n_thread, RHashtable **rhashtab, RHashtable **rrhashtab, uint32_t *sed_rec, uint16_t *sed_hit_times, 
	uint16_t *unused_bkt, opts *_opt, Aligner *aln) ;
static void 	*thread_worker(void *data);
extern 	uint8_t rev[];

#endif
