 /*
	Description: about the method what we used in select hit area 
	use 2048-len bucket



*/
#include "aligner.h"

#include "readfl.h"



#define GAPOPEN 2
#define GAPEXTENDED 1
#define MISMATCH 5
#define MATCH 1


#define SELECT_NUM 10


#define FORELEN  5

uint8_t rev[128]={
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,
	4,4,4,4,65,4,4,4,4,4,4,4,4,4,4,4,4,
	84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,4,
	4,4,4,65,4,4,4,4,4,4,4,4,4,4,4
};


//merge the five biggest area into a larger area first I will produce the five area.
//parameters: bkt2 
//algorithm detect if the next one is no more than one merge them and compare the choose the best one 

int compare_bkt2(const void *p,const void *q)
{
	const bkt2 *f = (bkt2 *)p;
	const bkt2 *t = (bkt2 *)q;
	return (int)(f->seq_num - t->seq_num);
}

int compare_sam(const void *p, const void *q, void *arg)
{
	int f = *(int *)p;
	int h = *(int *)q;
	Sam_Rec *s = (Sam_Rec *)arg;
 	return s[h].score - s[f].score;
}

Aligner::Aligner(opts *_opt)
{
	opt = _opt;
	sed_rec = new uint32_t[LEN_BASES];
	if(NULL == sed_rec) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);
	}
	
	sed_hit_times = new uint16_t[LEN_BASES];
	if (NULL == sed_hit_times) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);	
	}

	unused_bkt = new uint16_t[LEN_BASES];
	if (NULL == unused_bkt) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);	
	}

	RCRead = new char[opt->len_limit];
	if (NULL == RCRead) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);
	}


	//aligner = new Simpledp(len_sed,writecigar,cigarbuffer,LEN_LIMIT<<1,nextcigarBuffer,LEN_LIMIT<<1, &usedcigarsize);
	//if (NULL == aligner) {
	//	cout<<"Error when Applying for new space, now exit..."<<endl;
	//	exit(1);
	//}

}

Aligner::~Aligner()
{
	if (NULL != sed_rec)
		delete[] sed_rec;
	if (NULL != sed_hit_times)
		delete[] sed_hit_times;
	if (NULL != unused_bkt)
		delete[] unused_bkt;
	if (NULL != RCRead)
		delete[] RCRead;
}
//generate seed according to read 
//status: checked 
//problem make sure that length of bases is longer than LEN_BASES and len_sed is less than LEN_BASES 
uint32_t Aligner::gen_sed(char *bases, uint32_t len_bases)
{
	uint32_t pre;
	uint32_t count = 0;
	//initiate the first one and for the others will be obtained by the its previous value.
	//pay attention to limitation of the program there may be something wrong 
	// I think this phase is not safe enoguh since the length of bases may be not less than len_sed and the LEN_BASES - len_sed may be negative check this
	transfer(bases,&pre,opt->len_sed);
	//cout<<pre<<endl;
	uint32_t mask = 0xffffffff>>(32 - (opt->len_sed<<1));
	if((hashtab->pointer[pre]^hashtab->pointer[pre+1])!=0&&(hashtab->pointer[pre]&
		hashtab->pointer[pre+1])!=0) {
		sed_rec[0] = pre;
		++count;
	}

	for(uint32_t i=1;i<=len_bases - opt->len_sed;i++) {
		pre = ((pre<<2)&mask)|trans[bases[i + opt->len_sed - 1 ]];
		//cout<<pre<<endl;
		if((hashtab->pointer[pre]^hashtab->pointer[pre+1])!=0&&(hashtab->pointer[pre]&
		hashtab->pointer[pre+1])!=0) {
			sed_rec[count] = pre;
			++count;
		}
	}

	return count;
}

//status: checked 
//problem: none
int Aligner::RevComRead(char *read,int len_read)
{
	int i;
	for(i=0;i<len_read;i++)
		RCRead[i] = rev[read[len_read - 1 - i]];
	RCRead[i] = '\0';
	return 0;
}					


//status: checked
//problem: none
int compar(const void *p,const void *q)
{
	const uint32_t *t = (uint32_t *)p;
	const uint32_t *f = (uint32_t *)q;
	if(*t>*f)
		return 1;
	else
		if (*t<*f)
			return -1;
		else
			return 0;
}


void Aligner::sort_sed(uint32_t usedseed)
{
	qsort(sed_rec,usedseed,sizeof(uint32_t),compar);
}


void  Aligner::stat_sed(uint32_t usedseed)
{
	
	uint32_t j;
	for(uint32_t i=0;i<usedseed;) {
		for(j=1;i+j<usedseed&&sed_rec[i]==sed_rec[i+j];++j);//speed
			sed_hit_times[i] = j;
			unused_bkt[i] = hashtab->pointer[sed_rec[i]+1] - hashtab->pointer[sed_rec[i]];
			
			if (unused_bkt[i] >= opt->hit_limit)
				unused_bkt[i] = 0;
			i += j;		
	}
	
	
}

uint32_t fig_pos(uint32_t pos,uint32_t offset)
{
	uint32_t bkt_num = pos<<10;
	if (bkt_num <= offset) {
		return 0;
	}

	uint32_t left_start =  bkt_num - offset ; 
	return left_start;
}

void Aligner::proCans(char *read, uint32_t lenRead, bool isRC,  std::priority_queue<bkt2> &cansHeap) 
{				
	std::priority_queue<bkt> 	bkt_heap;
	std::priority_queue<bkt2> 	bkt_heap2;

	//uint32_t 		bkt_num;
	uint16_t 		bkt_hit_times, bkt_hit_times2;
	uint32_t 		bkt_seq_num_pre, bkt_seq_num_pre2;
	
	bkt 			bkt_temp;
	bkt2 			bkt_temp2;

	uint32_t 		pos_seq_num_pre;
	uint32_t 		array_seq_temp;


	uint32_t num_usedsed = gen_sed(read, lenRead);
	sort_sed(num_usedsed);
	stat_sed(num_usedsed);
	//cout<<"reverse one"<<endl;
	for(uint32_t i=0;i<num_usedsed;) {
		if (unused_bkt[i]) {
			bkt_temp.sequence_num_bkt = hashtab->seq_bkt[hashtab->pointer[sed_rec[i]] - 1];
			bkt_temp.array_seq = i;
			--unused_bkt[i];
			bkt_heap.push(bkt_temp);
		}
	 	i += sed_hit_times[i];
	}
	if (!bkt_heap.empty()) {
		pos_seq_num_pre = bkt_heap.top().sequence_num_bkt;
		array_seq_temp 	= bkt_heap.top().array_seq;
		
		//cout<<pos_seq_num_pre<<endl;

		bkt_hit_times 	= sed_hit_times[array_seq_temp];
		bkt_seq_num_pre = pos_seq_num_pre - 1;//zero


		bkt_hit_times2 	= bkt_hit_times;
		bkt_seq_num_pre2 = pos_seq_num_pre;
		
		bkt_heap.pop();
			//cout<<"pop:"<<bkt_seq_num_pre<<"  "<<array_seq_temp<<endl;

		if(0 != unused_bkt[array_seq_temp]) {
			bkt_temp.sequence_num_bkt = hashtab->seq_bkt[hashtab->pointer[sed_rec[array_seq_temp] + 1] - unused_bkt[array_seq_temp] -1];
			
			bkt_temp.array_seq = array_seq_temp;
			bkt_heap.push(bkt_temp);
				//cout<<"push:"<<bkt_temp.sequence_num_bkt<<"  "<<bkt_temp.array_seq<<endl;
			--unused_bkt[array_seq_temp];
		}
	}
	
	while(!bkt_heap.empty()) {		
			array_seq_temp = bkt_heap.top().array_seq;
		if(pos_seq_num_pre == bkt_heap.top().sequence_num_bkt) {//here may be use ^
			bkt_hit_times += sed_hit_times[array_seq_temp];
			//cout<<sed_hit_times[array_seq_temp]<<endl;
			bkt_hit_times2 += sed_hit_times[array_seq_temp];//pay attention to this place it is not bkt_hit_times2 = bkt_hit_times
			//cout<<bkt_hit_times<<"\t0"<<endl;
			
		} else {
			//cout<<bkt_seq_num_pre<<endl;
			if(pos_seq_num_pre + 1 == bkt_heap.top().sequence_num_bkt) {//here may be use ^
				if (~bkt_seq_num_pre) {
					bkt_temp2.seq_num 	= bkt_seq_num_pre;
					bkt_temp2.hit_times = bkt_hit_times;
					//cout<<bkt_hit_times<<"\t1"<<endl;
					bkt_temp2.isrc = isRC;
					bkt_heap2.push(bkt_temp2);	
				}
				++pos_seq_num_pre;
				bkt_hit_times2 += sed_hit_times[array_seq_temp];
				//cout<<sed_hit_times[array_seq_temp]<<"\t2"<<endl;
				bkt_hit_times = bkt_hit_times2;
				bkt_seq_num_pre = bkt_seq_num_pre2;

				bkt_hit_times2 = sed_hit_times[array_seq_temp];
				bkt_seq_num_pre2 = pos_seq_num_pre;
			} else {
				if (~bkt_seq_num_pre) {
					bkt_temp2.seq_num 	= bkt_seq_num_pre;
					bkt_temp2.hit_times = bkt_hit_times;
					//cout<<bkt_hit_times<<"\t2"<<endl;
					bkt_temp2.isrc = isRC;
					bkt_heap2.push(bkt_temp2);	
				}

				bkt_temp2.seq_num 	= bkt_seq_num_pre2;
				bkt_temp2.hit_times = bkt_hit_times2;
				//cout<<bkt_hit_times2<<"\t2"<<endl;
				bkt_temp2.isrc = isRC;
				bkt_heap2.push(bkt_temp2);	

				pos_seq_num_pre = bkt_heap.top().sequence_num_bkt;	

				bkt_hit_times = sed_hit_times[array_seq_temp];
				bkt_seq_num_pre = pos_seq_num_pre - 1;

				bkt_hit_times2 = bkt_hit_times;
				bkt_seq_num_pre2 = pos_seq_num_pre;
			}

		}
		//cout<<bkt_heap.top().sequence_num_bkt<<endl;
		bkt_heap.pop();
			//cout<<"pop:"<<bkt_seq_num_pre<<"  "<<array_seq_temp<<endl;
		if(0 != unused_bkt[array_seq_temp]) {
			bkt_temp.sequence_num_bkt = hashtab->seq_bkt[hashtab->pointer[sed_rec[array_seq_temp] + 1] - unused_bkt[array_seq_temp] -1];
			bkt_temp.array_seq = array_seq_temp;
			bkt_heap.push(bkt_temp);
			--unused_bkt[array_seq_temp];
		}
	}

	
	for(uint32_t i=0;i<opt->canN&&!bkt_heap2.empty();++i) {
		cansHeap.push(bkt_heap2.top());
		bkt_heap2.pop();
	}			
}
// for non sv 
int Aligner::conductAlign(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab
	,std::priority_queue <bkt2> &cansHeap) 
{
	
	uint32_t 		seq_num_temp;
	uint32_t 		diff = trunk->seq.l > LEN_BASES ? trunk->seq.l - LEN_BASES:0;
	uint32_t 		offset = diff>>1;
	int 			score;
	//uint32_t 		pos;
	uint32_t 		extract_length;
	bool 			rc;
	uint32_t 		left_start;
	
	bool 			hashinitiateRC = false;
	bool			hashinitiate = false;
	
	
	//int 			score;
	
	int 			flag = 0;
	int 			sign = 0;
	int 			usedArray = 0;

	//uint8_t 		posCount = 0;
	RHashtable 		*usedhash;
	char 			*_usedread;
	int 			countSam = 0;

	for (uint32_t i=0; i< opt->canN&&!cansHeap.empty();++i) {
		preserved[usedArray++] = cansHeap.top();

		Graphic 		gra;
		seq_num_temp = cansHeap.top().seq_num;
		rc = cansHeap.top().isrc;
		
		//if(~seq_num_temp) // what if there is not enough length
		//if (rc) useoffset = revoffset; else useoffset = offset;
		extract_length = diff + 2048;
		left_start = fig_pos(seq_num_temp,offset);
		if(left_start+extract_length>=len_genome) {
			extract_length = len_genome - left_start;
		} 
		
		if(rc) {
			if (!hashinitiateRC) {
				rrhashtab->buildRHash(RCRead, trunk->seq.l);		
				hashinitiateRC = true;
			}
			_usedread = RCRead;
			usedhash = rrhashtab;
					
		} else {
			if (!hashinitiate) {
				rhashtab->buildRHash(trunk->seq.s, trunk->seq.l);
				hashinitiate = true;
			}					
			_usedread = trunk->seq.s;
			usedhash = rhashtab;				//cout<<pos+left_start<<endl;	
		}			

		flag = gra.applyGraphic(usedhash, genome+left_start, extract_length, _usedread, trunk->seq.l,&score,opt->waitingLen, left_start,rc,Start_pos,ChrName,
				countChr, sams, countSam);

		if ( 1 == flag) {
			sign = 1;
			++countSam;
		} 
		
		cansHeap.pop();

	}
	
	if (1 != sign) {
		//recuculate
		for (int i=0;i<usedArray;++i) {
			Graphic 		gra;
			seq_num_temp = preserved[i].seq_num;
			rc = preserved[i].isrc;
			
			//if(~seq_num_temp) // what if there is not enough length
			//if (rc) useoffset = revoffset; else useoffset = offset;

			left_start = fig_pos(seq_num_temp,offset);
			if(left_start+extract_length>=len_genome) {
				extract_length = len_genome - left_start;
			} 
			
			if(rc) {

				if (!hashinitiateRC) {
					rrhashtab->buildRHash(RCRead, trunk->seq.l);	
					hashinitiateRC = true;	
				}
				_usedread = RCRead;
				usedhash = rrhashtab;
							
			} else {
				if (!hashinitiate) {
					rhashtab->buildRHash(trunk->seq.s, trunk->seq.l);
					hashinitiate = true;
				}
				_usedread = trunk->seq.s;
				usedhash = rhashtab;
						
			}

			flag = gra.applyGraphic(usedhash, genome+left_start, extract_length, _usedread, trunk->seq.l,&score,opt->waitingLen<<1, left_start,rc,Start_pos,ChrName,
					countChr, sams, countSam);
			if (flag == 1)  {
				++countSam;
				sign = 1;
			}		

		}

	}
	//output sam records
	if (sign == 1) {
		if (countSam >= 2) {
			int orders[countSam];
			for(int i=0;i<countSam;++i) orders[i] = i;
			qsort_r(orders,countSam,sizeof(int),compare_sam,sams);
			int quality;
			if (sams[orders[0]].score > 0) {
				quality = (int)(250.0 * 0.25 * (double)(sams[orders[0]].score - sams[orders[1]].score)/(double)(sams[orders[0]].score));
				if (quality>=60) quality = 60;
			} else quality = 0;
			if (sams[orders[0]].flag)
				 _usedread = RCRead;
			else
				_usedread = trunk->seq.s;

			cout<<trunk->name.s<<"\t"<<sams[orders[0]].flag<<"\t"<<ChrName[sams[orders[0]].chrIndex]<<"\t"<<sams[orders[0]].pos<<"\t"<<quality<<"\t"<<sams[orders[0]].headCigar<<sams[orders[0]].bodyCigar<<sams[orders[0]].tailCigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
			<<"0"<<"\t";
			cout<<_usedread<<"\t";
			cout<<trunk->qual.s<<"\t"<<"AS:i:"<<sams[orders[0]].score<<endl;

			for (int i=1;i<countSam;++i) {
				
				cout<<trunk->name.s<<"\t"<<sams[orders[i]].flag + 256<<"\t"<<ChrName[sams[orders[i]].chrIndex]<<"\t"<<sams[orders[i]].pos<<"\t"<<"0"<<"\t"<<sams[orders[i]].headCigar<<sams[orders[i]].bodyCigar<<sams[orders[i]].tailCigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
				<<"0"<<"\t";
				cout<<_usedread<<"\t";
				cout<<trunk->qual.s<<"\t"<<"AS:i:"<<sams[orders[i]].score<<endl;
			}


		} else {
			if (sams[0].flag)
				_usedread = RCRead;
			else
				_usedread = trunk->seq.s;
			cout<<trunk->name.s<<"\t"<<sams[0].flag<<"\t"<<ChrName[sams[0].chrIndex]<<"\t"<<sams[0].pos<<"\t"<<"60"<<"\t"<<sams[0].headCigar<<sams[0].bodyCigar<<sams[0].tailCigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
			<<"0"<<"\t";
			cout<<_usedread<<"\t";
			cout<<trunk->qual.s<<"\t"<<"AS:i:"<<sams[0].score<<endl;
		}
	}

	//calculate MAPQ
	

	while (!cansHeap.empty()) 
		cansHeap.pop();

	return sign;

} 

int Aligner::conductAlign(kseq_t *trunk, char *read, char *rcRead, int lenRead, RHashtable *rhashtab, RHashtable *rrhashtab
	,std::priority_queue <bkt2> &cansHeap, Sam_Rec *svsamsp) 
{
	
	uint32_t 		seq_num_temp;
	uint32_t 		diff = lenRead - LEN_BASES;
	uint32_t 		offset = diff>>1;
	int 			score;
	//uint32_t 		pos;
	uint32_t 		extract_length;
	bool 			rc;
	uint32_t 		left_start;
	
	bool 			hashinitiateRC = false;
	bool			hashinitiate = false;
	
	
	//int 			score;
	
	int 			flag = 0;
	int 			sign = 0;
	int 			usedArray = 0;

	//uint8_t 		posCount = 0;
	RHashtable 		*usedhash;
	char 			*_usedread;
	int 			countSam = 0;

	for (uint32_t i=0; i< opt->canN&&!cansHeap.empty();++i) {
		preserved[usedArray++] = cansHeap.top();

		Graphic 		gra;
		seq_num_temp = cansHeap.top().seq_num;
		rc = cansHeap.top().isrc;
		
		//if(~seq_num_temp) // what if there is not enough length
		//if (rc) useoffset = revoffset; else useoffset = offset;
		extract_length = diff + 2048;
		left_start = fig_pos(seq_num_temp,offset);
		if(left_start+extract_length>=len_genome) {
			extract_length = len_genome - left_start;
		} 
		
		if(rc) {
			if (!hashinitiateRC) {
				rrhashtab->buildRHash(rcRead, lenRead);		
				hashinitiateRC = true;
			}
			_usedread = rcRead;
			usedhash = rrhashtab;	
		} else {
			if (!hashinitiate) {
				rhashtab->buildRHash(read, lenRead);
				hashinitiate = true;
				
			}					
			_usedread = read;
			usedhash = rhashtab;				//cout<<pos+left_start<<endl;	
		}			

		flag = gra.applyGraphic(usedhash, genome+left_start, extract_length, _usedread, lenRead,&score,opt->waitingLen,left_start,rc,Start_pos,ChrName,
				countChr, svsamsp, countSam);

		if ( 1 == flag) {
			sign = 1;
			++countSam;
		} 
		
		cansHeap.pop();

	}
	
	while (!cansHeap.empty()) 
		cansHeap.pop();

	return countSam;

} 

int  Aligner::applyNonSV(kseq_t *trunk,  RHashtable *rhashtab, RHashtable *rrhashtab)
{
	
	std::priority_queue<bkt2> 	cansHeap;
	
	uint32_t diff = trunk->seq.l > LEN_BASES?trunk->seq.l - LEN_BASES:0;
	
	uint32_t offset = diff>>1;
	
	uint32_t canReadLen = trunk->seq.l > LEN_BASES ? LEN_BASES:trunk->seq.l;  

	char * useread = trunk->seq.s;
	bool isRC = false;
	//uint32_t extend = (uint32_t)(1.1*offset);//waiting to be tested;
	for (int j=0;j<=1;j++) {
		proCans(useread+offset, canReadLen, isRC, cansHeap);			
		useread = RCRead;
		isRC = true;
	}
	
	int sign = conductAlign(trunk,rhashtab,rrhashtab, cansHeap);
	sign = trunk->seq.l < LEN_BASES ? 1:sign; //if read length is less than len bases even though it was not mapped it won't conduct sv operation;
	return sign;
}

// void Aligner::applySV(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab)
// {
// 	//split into several pieces

// 	std::priority_queue<bkt2> 	cansHeap;
// 	uint32_t pNumber = trunk->seq.l/LEN_BASES;
// 	uint32_t leftLen = LEN_BASES + trunk->seq.l - pNumber*LEN_BASES;
// 	char *useread = trunk->seq.s;
// 	bool isRC = false;
// 	// produce the 
// 	for (uint32_t i = 0; i<pNumber-1;++i) {
// 		for (int j=0;j<=1;++j) {
// 			proCans(useread+i*LEN_BASES,isRC,cansHeap);
// 			useread = RCRead;
// 			isRC = true;
// 		}	
// 		conductAlign(trunk,trunk->seq.s+i*LEN_BASES, RCRead + i*LEN_BASES, LEN_BASES, rhashtab,rhashtab,cansHeap,i*LEN_BASES,trunk->seq.l-(i+1)*LEN_BASES);

// 	}
// 	//deal with the last piece	
// 	uint32_t offset = (leftLen - LEN_BASES)>>1;
// 	for (int j=0;j<=1;++j) {
// 		proCans(useread+(pNumber-1)*LEN_BASES + offset,isRC,cansHeap);
// 		useread = RCRead;
// 		isRC = true;
// 	}
	
// 	conductAlign(trunk,trunk->seq.s+(pNumber-1)*LEN_BASES,RCRead+(pNumber-1)*LEN_BASES,leftLen,rhashtab,rrhashtab,cansHeap,(pNumber-1)*LEN_BASES,0);

// }

//
int Aligner::Qualified(Sam_Rec *key,Sam_Rec *set, int start, int end, bool *issubstr)
{
	for (int i=start;i<end;++i) {
		if (!issubstr[i]) {
			if (key->chrIndex == set[i].chrIndex && key->flag == set[i].flag) {
				if (key->flag) {
					if ( set[i].ref_end <= key->ref_start  && set[i].ref_end + opt->waitingLen >= key->ref_start && set[i].read_end <= key->read_start && set[i].read_end + opt->waitingLen >= key->read_start)
						return i;
				} else {
					if (key->ref_end <= set[i].ref_start && key->ref_end + opt->waitingLen >= set[i].ref_start && key->read_end <= set[i].read_start && key->read_end + opt->waitingLen >= set[i].read_start) 
						return i;
				}
			}
		} 
	}
	return -1;
}
int 	transIntoDec(uint8_t *transtr,char *str,int length)
{
	for (int i=0;i<length;++i) {
		transtr[i] = seq_nt4_tablet[str[i]];
	}
	return 0;
}

int Aligner::connect(Sam_Rec *rec1, Sam_Rec *rec2, kseq_t *trunk) 
{
	// cout<<"sdfasdfds"<<endl;
	//rec1->ref_end = rec2->ref_end;
	//rec1->read_end = rec2->read_end;
	uint8_t readqry[1024];
	uint8_t refqry[1024];
	std::string transitionPart;
	char 	tempCigar[50];
	int 	bscore = 0;
	char *readStartP;
	char *refStartP;
	int read_len = rec2->read_start - rec1->read_end;
	//cout<<rec2->read_start<<"\t"<<rec1->read_end<<endl;
	//cout<<rec2->ref_start<<"\t"<<rec1->ref_end<<endl;
	int ref_len = rec2->ref_start - rec1->ref_end;
	// fprintf(stderr,"%hu\t%hu\n",rec2->read_start,rec1->read_end);
	// fprintf(stderr,"%hu\t%hu\n",rec2->ref_start,rec1->ref_end);
	// fprintf(stderr,"%d\t%d\n",ref_len,read_len);

	const 	char 	correspondTable[] = "MIDNSHP=X";
	if (read_len == 0 || ref_len == 0) {
		if (read_len != 0) {
			sprintf(tempCigar,"%dI",read_len);
			transitionPart.append(tempCigar);
			bscore = GAPOPEN + (read_len - 1)*GAPEXTENDED;
		} else {
			if (0 != ref_len) {
				sprintf(tempCigar,"%dD",ref_len);
				transitionPart.append(tempCigar);
				bscore = GAPOPEN + (ref_len - 1)*GAPEXTENDED;
			}
		}
	} else {

		if (rec1->flag) 
			readStartP = RCRead + rec1->read_end;
		else
			readStartP = trunk->seq.s + rec1->read_end;
		refStartP = genome + rec1->ref_start;
		transIntoDec(readqry,readStartP,read_len);
		transIntoDec(refqry,refStartP,ref_len);

		const uint8_t *readqry_ = readqry;
		const uint8_t *refqry_ = refqry;
		//prseq(readStartP,read_len,true);
		//prseq(refStartP,ref_len,true);
		uint32_t 	*cigar;
		int 		n_cigar = 0;
		int w = read_len > ref_len ? read_len : ref_len;
		
		bscore = ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,GAPOPEN,GAPEXTENDED,w,&n_cigar,&cigar);

		//fprintf(stderr,"%d %d %d %d %d %d\n",order[i-1],order[i],read_len,ref_len,score, n_cigar);
		for (int z=0;z<n_cigar;++z) {
			sprintf(tempCigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			//sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
			//++startPosCigar;
			transitionPart.append(tempCigar);
		}
		//cout<<temp<<"\t"<<score<<'\t'<<"4"<<endl;
		free(cigar);
	}
	

	rec1->bodyCigar += transitionPart + rec2->bodyCigar;
	rec1->tailCigar = rec2->tailCigar;
	rec1->bodyScore += bscore + rec2->bodyScore;
	rec1->tailScore = rec2->tailScore;
	rec1->score = rec1->headScore + rec1->bodyScore + rec1->tailScore;
	rec1->ref_end = rec2->ref_end;
	rec1->read_end = rec2->read_end;


}
//input parameters char * read, char * rcRead for multiple threads

int Aligner::produceSAM(int countbulks,int *sam4bulk, kseq_t *trunk, uint32_t *len)
{
	int quality = 0;
	bool besubstr[sam4bulk[countbulks]];
	int extend[sam4bulk[countbulks]];

	for (int i=0;i<sam4bulk[countbulks];++i) {
		besubstr[i] = false;
		extend[i] = 0;
	}
	
	for (int i=0;i<countbulks;++i) {
		for (int j=sam4bulk[i];j<sam4bulk[i+1];++j) {
			if (!besubstr[j]) {
				for (int m=i+1;m<countbulks;++m) {
					int q;
					if ((q=Qualified(svsams+j,svsams,sam4bulk[m],sam4bulk[m+1],besubstr))!= -1) {
						if (svsams[j].flag) {
							connect(svsams+q,svsams+j,trunk);
							extend[q] = extend[j] - 1;
							besubstr[j] = true;

						} else {
							connect(svsams+j, svsams+q, trunk);
							++extend[j];
							besubstr[q] = true;
						}
						
						//--sam4bulk[m];
					} else
						break;
				}
			}
			
			//used[j] = true;
		}
	}
	
	//statistic besubstr; if it is one
	int countOutput = 0; 
	for (int i=0;i<sam4bulk[countbulks];++i) {
		if (!besubstr[i]) ++countOutput;
	}

	int secbestscore,secbestInd,secbestbulk;
	int bestscore,bestInd,bestbulk;
	bestscore = 0x80000001;

	for (int i=0;i<countbulks;++i) {
		for (int j=sam4bulk[i];j<sam4bulk[i+1];++j)
		if (!besubstr[j]) {
			if (bestscore <= svsams[j].score) {
				secbestscore = bestscore;
				secbestInd = bestInd; 
				secbestbulk = bestbulk;
				bestscore = svsams[j].score;
				bestInd = j; 
				bestbulk = i;

			} else {
				if (secbestscore <= svsams[j].score) {
					secbestscore = svsams[j].score;
					secbestInd = j; 
					secbestbulk = i;
				}
			}
		}
	}
	if (countOutput == 1) quality = 60;
	else  {
			if (bestscore <= 0) quality = 0;
			else {
				int quality = (int)(250.0 * 0.25 * (double)(bestscore - secbestscore)/(double)(bestscore));
				if (quality>=60) quality = 60;
			}
	}
//output first:
	int rclip,lclip;
	char *usedread;
	if (svsams[bestInd].flag) {
		rclip = len[bestbulk+extend[bestInd]];
		lclip = trunk->seq.l - len[bestbulk+1];
		usedread = RCRead;
	} else {
		lclip = len[bestbulk];
		rclip = trunk->seq.l - len[bestbulk+1+extend[bestInd]];
		usedread = trunk->seq.s;
	}

	cout<<trunk->name.s<<"\t"<<svsams[bestInd].flag<<"\t"<<ChrName[svsams[bestInd].chrIndex]<<"\t"<<svsams[bestInd].pos<<"\t"<<quality<<"\t";

	if (lclip)  	cout<<lclip<<"S";
	cout<<svsams[bestInd].headCigar<<svsams[bestInd].bodyCigar<<svsams[bestInd].tailCigar;
	if (rclip) 	cout<<rclip<<"S";
	cout<<"\t"<<"*"<<"\t"<<"0"<<"\t"<<"0"<<"\t";
	
	cout<<trunk->qual.s<<"\t"<<"AS:i:"<<svsams[bestInd].score<<endl;

	if (countOutput!=1) {
		for (int i=0;i<countbulks;++i) {
			for (int j=sam4bulk[i];j<sam4bulk[i+1];++j) {
				if (!besubstr[j] && j != bestInd) {
					//output 
					if (svsams[j].flag) {
						rclip = len[i+extend[j]];
						lclip = trunk->seq.l - len[i+1];
						usedread = RCRead;
					} else {
						lclip = len[i];
						rclip = trunk->seq.l - len[i+1+extend[j]];
						usedread = trunk->seq.s;
					}

					cout<<trunk->name.s<<"\t"<<svsams[j].flag + 256<<"\t"<<ChrName[svsams[j].chrIndex]<<"\t"<<svsams[j].pos<<"\t"<<quality<<"\t";
					if (lclip)	cout<<lclip<<"S";
					cout<<svsams[j].headCigar<<svsams[j].bodyCigar<<svsams[j].tailCigar;
					if (rclip)	cout<<rclip<<"S";
					cout<<"\t"<<"*"<<"\t"<<"0"<<"\t"<<"0"<<"\t";
					cout<<usedread<<"\t";
					cout<<trunk->qual.s<<"\t"<<"AS:i:"<<svsams[j].score<<endl;
				}
			}
		}	
	}
}


void Aligner::applySV(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab)
{
	//split into several pieces
	//
	std::priority_queue<bkt2> 	cansHeap;
	uint32_t pNumber = trunk->seq.l/LEN_BASES;
	uint32_t leftLen = LEN_BASES + trunk->seq.l - pNumber*LEN_BASES;
	//if (leftLen) ++pNumber;
	// produce the
	int bkt_index[pNumber+1];
	uint32_t len[pNumber+1];
	bkt_index[0] = 0; 
	len[0] = 0;
	for (int i=0;i<pNumber -1;++i) len[i+1] = len[i] + LEN_BASES;
	len[pNumber] = len[pNumber-1] + leftLen;
	
	char *useread;
	bool isRC;
	int  movement ;

	for (uint32_t i = 0; i < pNumber - 1;++i) {
		useread = trunk->seq.s;
		isRC = false;
		movement = i*LEN_BASES;
		for (int j=0;j<=1;++j) {
			proCans(useread+movement,LEN_BASES, isRC,cansHeap);
			useread = RCRead;
			movement = trunk->seq.l - len[i+1]; 
			isRC = true;
		}	
		int samCounter = conductAlign(trunk,trunk->seq.s+i*LEN_BASES, RCRead + movement, LEN_BASES, rhashtab,rrhashtab,cansHeap,svsams+bkt_index[i]);
		//update svsams
		bkt_index[i+1] = bkt_index[i] + samCounter;
		for (int k=bkt_index[i];k<bkt_index[i+1];++k) {
			if (svsams[k].flag) {
				svsams[k].read_start += movement; // here remember the read length can only be less than 65535 otherwise the structure has to be changed.
				svsams[k].read_end += movement;

			} else {
				svsams[k].read_start += len[i];
				svsams[k].read_end += len[i];	
			}
			
		}
	}

	
	//deal with the last piece	
	//statistic first bucket hit times:
	//initiate hit_times all to zero
	uint32_t offset = ((leftLen - LEN_BASES)>>1);
	useread = trunk->seq.s;
	isRC = false;
	movement = (pNumber-1)*LEN_BASES;
	for (int j=0;j<=1;++j) {
		proCans(useread+ movement + offset, LEN_BASES, isRC,cansHeap);
		useread = RCRead;
		movement = 0;
		isRC = true;
	}

	//len[pNumber] = len[pNumber - 1] + leftLen;
	bkt_index[pNumber] = bkt_index[pNumber-1] + conductAlign(trunk,trunk->seq.s+(pNumber-1)*LEN_BASES,RCRead,leftLen,rhashtab,rrhashtab,cansHeap,svsams+bkt_index[pNumber-1]);
	for (int k=bkt_index[pNumber-1];k<bkt_index[pNumber];++k) {
		if (svsams[k].flag == 0) {
			svsams[k].read_start += len[pNumber-1];
			svsams[k].read_end += len[pNumber-1];
		} 
	}
	if (bkt_index[pNumber] <= 0) return ;
	produceSAM(pNumber, bkt_index, trunk, len);
}

// void Aligner::applySV(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab)
// {
// 	//split into several pieces
// 	//
// 	// first analysis reverse or not
// 	std::priority_queue<bkt2> 	cansHeap;
// 	uint32_t offset = ((trunk->seq.l - LEN_BASES)>>1);
// 	char 	*useread;
// 	bool 	isRC;
// 	int  	movement ;
// 	useread = trunk->seq.s;
// 	isRC = false;
// 	for (int j=0;j<=1;++j) {
// 		proCans(useread+offset,isRC,cansHeap);
// 		useread = RCRead;
// 		isRC = true;
// 	}
// 	if (cansHeap.top().isrc) {
// 		useread = RCRead;
// 		isRC = true;
// 	} else {
// 		useread = trunk->seq.s;
// 		isRC = false;
// 	}
// 	while (!cansHeap.empty()) cansHeap.pop();

	
// 	uint32_t pNumber = trunk->seq.l/LEN_BASES;
// 	uint32_t leftLen = LEN_BASES + trunk->seq.l - pNumber*LEN_BASES;
// 	//if (leftLen) ++pNumber;
// 	// produce the
// 	int bkt_index[pNumber+1];
// 	uint32_t len[pNumber+1];
// 	bkt_index[0] = 0; 
// 	len[0] = 0;

// 	for (int i=0;i<pNumber -1;++i) len[i+1] = len[i] + LEN_BASES;
// 	len[pNumber] = len[pNumber-1] + leftLen;
	
	
// 	//
// 	for (uint32_t i = 0; i < pNumber - 1;++i) {
// 		movement = i*LEN_BASES;
// 		proCans(useread+movement,isRC,cansHeap);
// 		int samCounter = conductAlign(trunk,trunk->seq.s+i*LEN_BASES, RCRead + movement, LEN_BASES, rhashtab,rrhashtab,cansHeap,svsams+bkt_index[i]);
// 		//update svsams
// 		bkt_index[i+1] = bkt_index[i] + samCounter;
// 		for (int k=bkt_index[i];k<bkt_index[i+1];++k) {
// 			svsams[k].read_start += len[i];
// 			svsams[k].read_end += len[i];	
// 		}
// 	}
	
	
// 	//deal with the last piece	
// 	//statistic first bucket hit times:
// 	//initiate hit_times all to zero
// 	offset = ((leftLen - LEN_BASES)>>1);

// 	movement = (pNumber-1)*LEN_BASES;

// 	proCans(useread+ movement + offset,isRC,cansHeap);
		
// 	//len[pNumber] = len[pNumber - 1] + leftLen;
// 	bkt_index[pNumber] = bkt_index[pNumber-1] + conductAlign(trunk,trunk->seq.s+(pNumber-1)*LEN_BASES,RCRead+movement,leftLen,rhashtab,rrhashtab,cansHeap,svsams+bkt_index[pNumber-1]);
	

// 	for (int k=bkt_index[pNumber-1];k<bkt_index[pNumber];++k) {
// 		//if (svsams[k].flag == 0) {
// 		svsams[k].read_start += len[pNumber-1];
// 		svsams[k].read_end += len[pNumber-1];
		
// 	}
// 	if (bkt_index[pNumber] <= 0) return ;
// 	produceSAM(pNumber, bkt_index, trunk, len);
	
// }
void Aligner::Runtask()
{
	//read datas
	len_genome = 0;
	read_file rd;
	
	ChrName = new char *[LEN];
	for (int i=0;i<100;++i) { ChrName[i] = new char[LEN];}
	//int 	 countChr;
	genome = rd.read_ref(opt->refpath,&len_genome,Start_pos,ChrName,&countChr);
	
	if (NULL == genome) {
		
		fprintf(stderr,"Fail to load Genome reference, now exit");
		exit(1);
	} 

	//output header
	cout<<"@HD\tVN:"<<"v.15.01"<<endl;
	for (int i=1;i<=countChr;++i) {cout<<"@SQ\tSN:"<<ChrName[i]<<"\tLN:"<<Start_pos[i]-Start_pos[i-1]<<endl;}
	cout<<"@PG\tID:rHAT-mapper\tVN:v.15.01\tCL:";
	for (int i=0;i<opt->argc;++i) {cout<<opt->argv[i]<<" ";}
	cout<<endl;	

	Hash hashh;
	hashtab = hashh.load_hashfile(opt->hashdir,len_genome,opt->len_sed);

	if ( NULL == hashtab) {
		fprintf(stderr,"Fail to load Hash Table, now exit");
		exit(1);
	} 
	
	gzFile fp;
	kseq_t *trunk;
	
	fp = gzopen(opt->readpath, "r");
	trunk = kseq_init(fp);
	
	
	RHashtable *rhashtab = new RHashtable(opt->rh_seed_len,FORELEN,opt->len_limit);
	if ( NULL==rhashtab ) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

	//Graphic gra = new Graphic();

	//if (NULL == gra) { printf("can't allocate space to graphic now exit(1)!"); exit(1);}
	RHashtable *rrhashtab = new RHashtable(opt->rh_seed_len,FORELEN,opt->len_limit);
	if ( NULL== rrhashtab) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

	//lv = new LandauVishkin(13);

	//if (NULL == lv)  { printf("can't allocate space now exit"); exit(1);}
	//apply spaces for sam records
	
	//sams = (Sam_Rec *)malloc(sizeof(Sam_Rec)*opt->canN);
	sams = new Sam_Rec[opt->canN];
	if (NULL == sams) {fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}
	svsams = new Sam_Rec[opt->canN*((opt->len_limit/LEN_BASES)+1)];
	if (NULL == svsams) {fprintf(stderr, "Failed when applying for new space! now exit");exit(1);} 
	while(kseq_read(trunk) >= 0) {
		RevComRead(trunk->seq.s,trunk->seq.l);
		int sign = applyNonSV(trunk,rhashtab,rrhashtab);
		if (sign!=1) {
			applySV(trunk,rhashtab,rrhashtab);
		}
	}
	//cout<<"total number: "<<real_num<<endl;

	//delete gra;
	
	if (NULL != genome)
		delete[] genome;

	if(NULL != hashtab->pointer)
		delete[] hashtab->pointer;

	/*
	for (uint32_t i=0; i<=len_genome - len_sed; ++i) {
	    if( NULL != hashtab->seq_bkt[i])
	        delete[] hashtab->seq_bkt[i]; 
	}
	*/
	if (NULL!=hashtab->seq_bkt)
	    delete[] hashtab->seq_bkt;
	


	if(NULL != hashtab)
		delete hashtab;
	
	
	if (NULL != ChrName) {
		for (int i=0;i<100;++i) { delete ChrName[i]; }
		delete ChrName;
	}

	delete rhashtab;
	delete rrhashtab;
	kseq_destroy(trunk);
	gzclose(fp);
	
	
}
