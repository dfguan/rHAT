#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
using namespace std;


#include "hash.h"
#include "ksw.h"

typedef struct sam_rec{
	uint8_t 	chrIndex;
	uint32_t 	pos;
	//char 		_cigar[LEN_LIMIT<<1];
	//string		headCigar;
	//string 		bodyCigar;
	//string 		tailCigar;
	string		cigar;
	int 		score;

	//int 		headScore;
	//int 		bodyScore;
	//int 		tailScore;
	//int 		score;
	uint16_t 	flag;
	int16_t		MAQ;
	struct sam_rec &	operator=(const struct sam_rec &r)  
	{ 
		chrIndex = r.chrIndex; 
		pos = r.pos;
		cigar = r.cigar;
		score = r.score;
		flag = r.flag;
		MAQ = r.MAQ;
		return *this;
	}
}Sam_Rec;


typedef struct svsam_rec {
	uint8_t 	chrIndex;
	uint32_t 	pos;
	//char 		_cigar[LEN_LIMIT<<1];
	string		cigar;
	int 		score;
	uint16_t 	flag;
	
	uint32_t 	ref_end;
	uint16_t 	read_end;
	uint16_t	read_start;
	uint32_t 	ref_start;
	uint16_t	lclip;
	uint16_t 	rclip;
	int16_t 	MAQ;
	struct svsam_rec &	operator=(const struct svsam_rec &r)  
	{ 
		chrIndex = r.chrIndex; 
		pos = r.pos;
		// headCigar = r.headCigar;
		// bodyCigar = r.bodyCigar;
		// tailCigar = r.tailCigar;
		cigar = r.cigar;
		
		score = r.score;
		flag = r.flag;
		ref_start = r.ref_start;
		rclip = r.rclip;
		lclip = r.lclip;
		MAQ = r.MAQ;
		return *this;
	}
}SvSam_Rec;
// typedef struct {
// 	uint8_t chrIndex;
// 	uint32_t pos;
// 	char 	_cigar[LEN_BASES<<2];
// 	int 	score;
// 	uint16_t flag;
// 	uint32_t refEnd;
// }SvSam_Rec;

class RHashtable {
public:
	uint16_t 		*p2leftSeq;
	uint16_t 		*left_seq;
	uint16_t 		*p2seqNum;
	uint16_t  		*seq_num;
	uint32_t 		len_sed;
	uint32_t 		forelength;
	//char 				*readseq;
private:	
	uint32_t 		*kmer_value;
	//uint16_t		*seq_number;
	//uint16_t 		*order;	
	uint32_t 		limitOfp2leftSeq;
	

public:
	RHashtable(uint32_t kmer,uint32_t forelen, uint32_t len_limit);
	//RHashtable(uint32_t kmer,uint32_t forelen, uint32_t len_limit, char *re);//for debug
	~RHashtable();
	void buildRHash(char *seq, uint32_t lenSeq);
};



class vertex{
public:
    uint16_t		read_seq;
    uint16_t 		ref_seq;
    uint16_t 		len;
    bool 	operator<(const vertex &r) const { return read_seq < r.read_seq; }
};

class ASeed{
public:
    uint16_t		read_seq;
   // uint16_t 		ref_seq;
    uint16_t 		left_time;
    bool 	operator<(const ASeed &r) const { return read_seq < r.read_seq; }
    ASeed &	operator=(const ASeed &r)  { read_seq = r.read_seq; left_time = r.left_time; return *this;}
    //bool 	operator == (const ASeed &r) const { return ((read_seq == r.read_seq)&&(ref_seq == r.ref_seq));}
};



class Graphic{
private:
	vector<vertex> 	node;
    vector<ASeed> 	livingSeed;
    uint8_t			readqry[1024];
    uint8_t 		refqry[1024];
   	char 			revreadqry[1024];
   	char 			revrefqry[1024];
   	char 			*ref_s;
   	char 			*ref_t;
    
public:
			Graphic(char *_ref_s, char *_ref_e) ;
	int 	applyGraphic(RHashtable *rhashtab, char *ref, uint32_t lenRef, char *read, uint32_t lenRead,int *score, uint32_t waitingLen, uint32_t left_start,
		bool rc, uint32_t *startPos, char **chrName, int countChr,Sam_Rec *sam, int countSam, int8_t *mat, int gapo, int gape);
	int 	applyGraphic(RHashtable *rhashtab, char *ref, uint32_t lenRef, char *read, uint32_t lenRead,int *score, uint32_t waitingLen, uint32_t left_start,
		bool rc, uint32_t *startPos, char **chrName, int countChr,SvSam_Rec *sam, int countSam, int8_t *mat, int gapo, int gape);
private:
	int 	transIntoDec(uint8_t *transstr,char *str,int length);
	void 	buildCounter(char *seq, uint32_t len_seq, RHashtable *rhashtab,uint16_t *seq_counter, uint16_t *p2startPos);
	//void	buildCounter(char *seq, uint32_t len_seq, RHashtable *rhashtab,uint16_t *seq_counter, uint16_t *p2startPos, char *readdebug);
	int 	createVertex(uint16_t *seq_n, char *ref, uint16_t *seq_counter, uint16_t *p2startPos, uint32_t lenRef, uint32_t offset_ref, uint32_t kmer, char *read,uint32_t lenRead, uint32_t offset_read, uint32_t vertex_lim);
	int 	CalEditDistancewithCigar(int *order, int order_len, char *read, uint32_t readlen, char *ref, uint32_t reflen, uint32_t left_start,
			bool rc, uint32_t *startPos, char **chrName, int countChr,Sam_Rec *sam, int countSam, int8_t *mat, int gapo, int gape);
	int 	CalEditDistancewithCigar(int *order, int order_len, char *read, uint32_t readlen, char *ref, uint32_t reflen, uint32_t left_start,
			bool rc, uint32_t *startPos, char **chrName, int countChr,SvSam_Rec *sam, int countSam, int8_t *mat, int gapo, int gape);
	void 	revstr(uint8_t *revstring, char *string, int len);
	void 	dealCigar(char *cigar, char *headbuf, int headbuflen);
	int 	findPos(uint32_t lenRef, uint32_t lenRead, uint32_t waitingLen,bool type, vertex *vnode);
	void 	buildGraphic(RHashtable *rhashtab, char *seq, uint32_t lenSeq, char *read, uint32_t lenRead);
	int 	dealGraph(uint32_t lenRef, uint32_t lenRead, char *read, char *ref, int *score, uint32_t waitingLen, uint32_t left_start,
			bool rc, uint32_t *startPos, char **chrName, int countChr,Sam_Rec *sam, int countSam, int8_t *matrix, int gapo, int gape);
	int 	dealGraph(uint32_t lenRef, uint32_t lenRead, char *read, char *ref, int *score, uint32_t waitingLen, uint32_t left_start,
			bool rc, uint32_t *startPos, char **chrName, int countChr,SvSam_Rec *sam, int countSam, int8_t *matrix, int gapo, int gape);
	void 	createLimVertex(uint16_t *seq_n,char *ref, uint32_t lenRef, char *read, uint32_t lenRead, uint16_t *seq_counter,uint16_t *p2startPos, uint32_t kmer, vertex head, vertex tail);
	//void 	buildEdge(Graphic *gra);
};
extern const uint8_t seq_nt4_tablet[];
