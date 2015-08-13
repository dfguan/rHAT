 /*
	Description: about the method what we used in select hit area use 2048-len bucket 
*/

#include "aligner.h" 

#include "readfl.h" 

#include <pthread.h> 

#define SELECT_NUM 10

#define N_NEEDED 5000

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

int read_seq;
pthread_rwlock_t rwlock;

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
	sed_rec = new uint32_t[opt->thread * LEN_BASES];
	if(NULL == sed_rec) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);
	}

	sed_hit_times = new uint16_t[opt->thread * LEN_BASES];
	if (NULL == sed_hit_times) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);
	}

	unused_bkt = new uint16_t[opt->thread * LEN_BASES];
	if (NULL == unused_bkt) {
		fprintf(stderr, "Failed when applying for new space! now exit");
		exit(1);
	}

	// RCRead = new char[opt->thread * opt->len_limit];
	// if (NULL == RCRead) {
	// 	fprintf(stderr, "Failed when applying for new space! now exit");
	// 	exit(1);
	// }

	//initiate ksw parameters
	int k;
	for (int i = k = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
			mat[k++] = i == j? opt->match : -opt->mismatch;
		mat[k++] = 0; // ambiguous base
	}
	for (int j = 0; j < 5; ++j) mat[k++] = 0;

	//

	//lv = new LandauVishkin(13);

	//if (NULL == lv)  { printf("can't allocate space now exit"); exit(1);}
	//apply spaces for sam records

	//sams = (Sam_Rec *)malloc(sizeof(Sam_Rec)*opt->canN);
	//sams = new Sam_Rec[opt->canN * opt->thread];
	//if (NULL == sams) {fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

	//sv_interval = opt->canN * ((opt->len_limit/LEN_BASES)+1) ;
	//svsams = new Sam_Rec[ sv_interval * opt->thread];
	//if (NULL == svsams) {fprintf(stderr, "Failed when applying for new space! now exit");exit(1);}
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
	// if (NULL != RCRead)
	// 	delete[] RCRead;
}
//generate seed according to read
//status: checked
//problem make sure that length of bases is longer than LEN_BASES and len_sed is less than LEN_BASES
uint32_t Aligner::gen_sed(char *bases, uint32_t len_bases, uint32_t *_sed_rec)
{
	uint32_t pre;
	uint32_t count = 0;
	//uint32_t offset = groupNum * LEN_BASES;
	//initiate the first one and for the others will be obtained by the its previous value.
	//pay attention to limitation of the program there may be something wrong
	// I think this phase is not safe enoguh since the length of bases may be not less than len_sed and the LEN_BASES - len_sed may be negative check this
	transfer(bases,&pre,opt->len_sed);
	//cout<<pre<<endl;
	uint32_t mask = 0xffffffff>>(32 - (opt->len_sed<<1));
	if((hashtab->pointer[pre]^hashtab->pointer[pre+1])!=0&&(hashtab->pointer[pre]&
		hashtab->pointer[pre+1])!=0) {
		_sed_rec[0] = pre;
		++count;
	}

	for(uint32_t i=1;i<=len_bases - opt->len_sed;i++) {
		pre = ((pre<<2)&mask)|trans[bases[i + opt->len_sed - 1 ]];
		//cout<<pre<<endl;
		if((hashtab->pointer[pre]^hashtab->pointer[pre+1])!=0&&(hashtab->pointer[pre]&
		hashtab->pointer[pre+1])!=0) {
			_sed_rec[count] = pre;
			++count;
		}
	}

	return count;
}

//status: checked
//problem: none
int Aligner::revComRead(char *read, char *rcRead, int len_read)
{
	int i;
	//uint32_t offset = groupNum * opt->len_limit;
	for(i=0;i<len_read;i++)
		rcRead[i] = rev[read[len_read - 1 - i]];
	rcRead[i] = '\0';
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


void Aligner::sort_sed(uint32_t *_sed_rec,uint32_t usedseed)
{
	//uint32_t offset = groupNum * LEN_BASES;
	qsort(_sed_rec,usedseed,sizeof(uint32_t),compar);
}


void  Aligner::stat_sed(uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt, uint32_t usedseed)
{

	uint32_t j;
	//uint32_t offset = groupNum * LEN_BASES;

	for(uint32_t i=0;i<usedseed;) {
		for(j=1;i+j<usedseed&&_sed_rec[i]==_sed_rec[i+j];++j);//speed
			_sed_hit_times[i] = j;
			_unused_bkt[i] = hashtab->pointer[_sed_rec[i]+1] - hashtab->pointer[_sed_rec[i]];

			if (_unused_bkt[i] >= opt->hit_limit)
				_unused_bkt[i] = 0;
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

void Aligner::proCans(char *read, uint32_t lenRead, bool isRC, std::priority_queue<bkt2> &cansHeap,
	uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt)
{
	std::priority_queue<bkt> 	bkt_heap;
	std::priority_queue<bkt2> 	bkt_heap2;
	//uint32_t 		offset = groupNum * LEN_BASES;
	//uint32_t 		bkt_num;
	uint16_t 		bkt_hit_times, bkt_hit_times2;
	uint32_t 		bkt_seq_num_pre, bkt_seq_num_pre2;

	bkt 			bkt_temp;
	bkt2 			bkt_temp2;

	uint32_t 		pos_seq_num_pre;
	uint32_t 		array_seq_temp;


	uint32_t num_usedsed = gen_sed(read, lenRead, _sed_rec);
	sort_sed(_sed_rec, num_usedsed);
	stat_sed(_sed_rec,_sed_hit_times,_unused_bkt,num_usedsed);
	//cout<<"reverse one"<<endl;
	for(uint32_t i=0;i<num_usedsed;) {
		if (_unused_bkt[i]) {
			bkt_temp.sequence_num_bkt = hashtab->seq_bkt[hashtab->pointer[_sed_rec[i]] - 1];
			bkt_temp.array_seq = i;
			--_unused_bkt[i];
			bkt_heap.push(bkt_temp);
		}
	 	i += _sed_hit_times[i];
	}
	if (!bkt_heap.empty()) {
		pos_seq_num_pre = bkt_heap.top().sequence_num_bkt;
		array_seq_temp 	= bkt_heap.top().array_seq;

		//cout<<pos_seq_num_pre<<endl;

		bkt_hit_times 	= _sed_hit_times[array_seq_temp];
		bkt_seq_num_pre = pos_seq_num_pre - 1;//zero


		bkt_hit_times2 	= bkt_hit_times;
		bkt_seq_num_pre2 = pos_seq_num_pre;

		bkt_heap.pop();
			//cout<<"pop:"<<bkt_seq_num_pre<<"  "<<array_seq_temp<<endl;

		if(0 != _unused_bkt[array_seq_temp]) {
			bkt_temp.sequence_num_bkt = hashtab->seq_bkt[hashtab->pointer[_sed_rec[array_seq_temp] + 1] - _unused_bkt[array_seq_temp] -1];

			bkt_temp.array_seq = array_seq_temp;
			bkt_heap.push(bkt_temp);
				//cout<<"push:"<<bkt_temp.sequence_num_bkt<<"  "<<bkt_temp.array_seq<<endl;
			--_unused_bkt[array_seq_temp];
		}
	}

	while(!bkt_heap.empty()) {
			array_seq_temp = bkt_heap.top().array_seq;
		if(pos_seq_num_pre == bkt_heap.top().sequence_num_bkt) {//here may be use ^
			bkt_hit_times += _sed_hit_times[array_seq_temp];
			//cout<<sed_hit_times[array_seq_temp]<<endl;
			bkt_hit_times2 += _sed_hit_times[array_seq_temp];//pay attention to this place it is not bkt_hit_times2 = bkt_hit_times
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
				bkt_hit_times2 += _sed_hit_times[array_seq_temp];
				//cout<<sed_hit_times[array_seq_temp]<<"\t2"<<endl;
				bkt_hit_times = bkt_hit_times2;
				bkt_seq_num_pre = bkt_seq_num_pre2;

				bkt_hit_times2 = _sed_hit_times[array_seq_temp];
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

				bkt_hit_times = _sed_hit_times[array_seq_temp];
				bkt_seq_num_pre = pos_seq_num_pre - 1;

				bkt_hit_times2 = bkt_hit_times;
				bkt_seq_num_pre2 = pos_seq_num_pre;
			}

		}
		//cout<<bkt_heap.top().sequence_num_bkt<<endl;
		bkt_heap.pop();
			//cout<<"pop:"<<bkt_seq_num_pre<<"  "<<array_seq_temp<<endl;
		if(0 != _unused_bkt[array_seq_temp]) {
			bkt_temp.sequence_num_bkt = hashtab->seq_bkt[hashtab->pointer[_sed_rec[array_seq_temp] + 1] - _unused_bkt[array_seq_temp] -1];
			bkt_temp.array_seq = array_seq_temp;
			bkt_heap.push(bkt_temp);
			--_unused_bkt[array_seq_temp];
		}
	}


	for(uint32_t i=0;i<opt->canN&&!bkt_heap2.empty();++i) {
		cansHeap.push(bkt_heap2.top());
		bkt_heap2.pop();
	}
}
// for non sv
int Aligner::conductAlign(kseq_t *trunk,std::priority_queue <bkt2> &cansHeap, RHashtable *_rhashtable, RHashtable *_rrhashtable,
	Sam_Rec *_sams)
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
	//uint32_t 		offset = groupNum * opt->len_limit;

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

		Graphic 		gra(genome, genome_e);
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
				_rrhashtable->buildRHash(trunk->seq.rs, trunk->seq.l);
				hashinitiateRC = true;
			}
			_usedread = trunk->seq.rs;
			usedhash = _rrhashtable;

		} else {
			if (!hashinitiate) {
				_rhashtable->buildRHash(trunk->seq.s, trunk->seq.l);
				hashinitiate = true;
			}
			_usedread = trunk->seq.s;
			usedhash = _rhashtable;				//cout<<pos+left_start<<endl;
		}

		flag = gra.applyGraphic(usedhash, genome+left_start, extract_length, _usedread, trunk->seq.l,&score,opt->waitingLen, left_start,rc,Start_pos,ChrName,
				countChr, _sams, countSam, mat, opt->gapopen, opt->gapextend);

		if ( 1 == flag) {
			sign = 1;
			++countSam;
		}

		cansHeap.pop();

	}
	bool extendWaitingLen = false;

	if (1 != sign) {
		extendWaitingLen = true;
		//recuculate extending waiting length
		for (int i=0;i<usedArray;++i) {
			Graphic 		gra(genome, genome_e);
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
					_rrhashtable->buildRHash(trunk->seq.rs, trunk->seq.l);
					hashinitiateRC = true;
				}
				_usedread = trunk->seq.rs;
				usedhash = _rrhashtable;

			} else {
				if (!hashinitiate) {
					_rhashtable->buildRHash(trunk->seq.s, trunk->seq.l);
					hashinitiate = true;
				}
				_usedread = trunk->seq.s;
				usedhash = _rhashtable;

			}

			flag = gra.applyGraphic(usedhash, genome+left_start, extract_length, _usedread, trunk->seq.l,&score,opt->waitingLen<<1, left_start,rc,Start_pos,ChrName,
					countChr, _sams, countSam, mat, opt->gapopen, opt->gapextend);
			if (flag == 1)  {
				++countSam;
				sign = 1;
			}

		}

	}
	if (sign == 1) {
		int quality;
		if (countSam >= 2) {
			int orders[countSam];
			if (extendWaitingLen) {
				quality = 0;
			} else {
				for(int i=0;i<countSam;++i) orders[i] = i;
				qsort_r(orders,countSam,sizeof(int),compare_sam,_sams);
				if (_sams[orders[0]].score > 0) {
					quality = (int)(250.0 * 0.25 * (double)(_sams[orders[0]].score - _sams[orders[1]].score)/(double)(_sams[orders[0]].score));
					if (quality>=60) quality = 60;
				} else quality = 0;
				// if (sams[orders[0]].flag)
				// 	 _usedread = RCRead;
				// else
				// 	_usedread = trunk->seq.s;
				_sams[orders[0]].MAQ = quality;
				// cout<<trunk->name.s<<"\t"<<sams[orders[0]].flag<<"\t"<<ChrName[sams[orders[0]].chrIndex]<<"\t"<<sams[orders[0]].pos<<"\t"<<quality<<"\t"<<sams[orders[0]].headCigar<<sams[orders[0]].bodyCigar<<sams[orders[0]].tailCigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
				// <<"0"<<"\t";
				// cout<<_usedread<<"\t";
				// cout<<usedqual<<"\t"<<"AS:i:"<<sams[orders[0]].score<<endl;

				for (int i=1;i<countSam;++i) {
					_sams[orders[i]].flag += 256;
					_sams[orders[i]].MAQ = 0;

					// cout<<trunk->name.s<<"\t"<<sams[orders[i]].flag + 256<<"\t"<<ChrName[sams[orders[i]].chrIndex]<<"\t"<<sams[orders[i]].pos<<"\t"<<"0"<<"\t"<<sams[orders[i]].headCigar<<sams[orders[i]].bodyCigar<<sams[orders[i]].tailCigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
					// <<"0"<<"\t";
					// cout<<_usedread<<"\t";
					// cout<<usedqual<<"\t"<<"AS:i:"<<sams[orders[i]].score<<endl;
				}
				//may be if shoud be added
				Sam_Rec temp = _sams[orders[0]];
				_sams[orders[0]] = _sams[0];
				_sams[0] = temp;
			}
		}  else  {

		 _sams[0].MAQ = 60;

		}
	}
	while (!cansHeap.empty())
		cansHeap.pop();

	return countSam;

}

int Aligner::outputSam(kseq_t *_seqs, Sam_Rec *_sams, SvSam_Rec **_svsams, uint16_t *_sam_details, int _n_seqs)
{
	char *probQual = "*";
	char *usedqual;
	char *_usedread;
	kseq_t *trunk;

	for (int i=0; i<_n_seqs; ++i) {
			if (_sam_details[i]) {
				int countSam = _sam_details[i] >> 1;
				trunk = _seqs + i;
				if (_sam_details[i]&1) {
					Sam_Rec *nonSv = _sams + i * opt->canN;
					for (int j=0; j< countSam; ++j) {
						if (trunk->qual.s == NULL)// this might happend so qual is wrong?
							usedqual = probQual;
						else
							usedqual = trunk->qual.s;

						if (nonSv[j].flag)
							_usedread = trunk->seq.rs;
						else
							_usedread = trunk->seq.s;

						cout<<trunk->name.s<<"\t"<<nonSv[j].flag<<"\t"<<ChrName[nonSv[j].chrIndex]<<"\t"<<nonSv[j].pos<<"\t"<<nonSv[j].MAQ<<"\t"<<nonSv[j].cigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
						<<"0"<<"\t";
						cout<<_usedread<<"\t";
						cout<<usedqual<<"\t"<<"AS:i:"<<nonSv[j].score<<endl;

					}
				} else {
					//output svsam
					SvSam_Rec *sv = _svsams[i];
                    //cout<<"splitup"<<endl;
					for (int j=0; j < countSam; ++j) {

						if (trunk->qual.s == NULL)// this might happend so qual is wrong?
							usedqual = probQual;
						else
							usedqual = trunk->qual.s;
						if (sv[j].flag) 		_usedread = trunk->seq.rs;
						else 					_usedread = trunk->seq.s;

						cout<<trunk->name.s<<"\t"<<sv[j].flag<<"\t"<<ChrName[sv[j].chrIndex]<<"\t"<<sv[j].pos<<"\t"<<sv[j].MAQ<<"\t";
						if (sv[j].lclip)	cout<<sv[j].lclip<<"S";
						cout<<sv[j].cigar;
						if (sv[j].rclip)	cout<<sv[j].rclip<<"S";
						cout<<"\t"<<"*"<<"\t"<<"0"<<"\t"<<"0"<<"\t";
						cout<<_usedread<<"\t";
						cout<<usedqual<<"\t"<<"AS:i:"<<sv[j].score<<endl;


					}

				}

			}
		}
		return 0;

}

int Aligner::conductAlign(kseq_t *trunk, char *read, char *rcRead, int lenRead, RHashtable *_rhashtab, RHashtable *_rrhashtab
	,std::priority_queue <bkt2> &cansHeap, SvSam_Rec *_svsamsp)
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

		Graphic 		gra(genome, genome_e);
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
				_rrhashtab->buildRHash(rcRead, lenRead);
				hashinitiateRC = true;
			}
			_usedread = rcRead;
			usedhash = _rrhashtab;
		} else {
			if (!hashinitiate) {
				_rhashtab->buildRHash(read, lenRead);
				hashinitiate = true;

			}
			_usedread = read;
			usedhash = _rhashtab;				//cout<<pos+left_start<<endl;
		}

		flag = gra.applyGraphic(usedhash, genome+left_start, extract_length, _usedread, lenRead,&score,opt->waitingLen,left_start,rc,Start_pos,ChrName,
				countChr, _svsamsp, countSam, mat, opt->gapopen, opt->gapextend);

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



int Aligner::qualified(SvSam_Rec *key,SvSam_Rec *set, int start, int end)
{
	for (int i=start;i<end;++i) {
		if (~set[i].MAQ) {
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
int Aligner::mergeBoundaryCigar(string & cigar1, string & cigar2, uint32_t *cigar, int n_cigar, const char *correspondTable) 
{
	uint32_t len_cigar1 = cigar1.size();
	uint32_t len_cigar2 = cigar2.size();
	
	uint32_t digit_end_cigar2,digit_start_cigar1;

	for (digit_end_cigar2 = 0; digit_end_cigar2 < len_cigar2; ++digit_end_cigar2) {
		if (isalpha(cigar2[digit_end_cigar2]))
			break;
	}

	for (digit_start_cigar1 = len_cigar1 -2; digit_start_cigar1 >0;--digit_start_cigar1) {
		if (isalpha(cigar1[digit_start_cigar1])) 
			break;
	}

	uint32_t len_digit = len_cigar1 - digit_start_cigar1 - 1; 
	uint32_t countAlpha_cigar1 = strtoul(cigar1.substr(digit_start_cigar1 + 1, len_digit).c_str(),0,10);
	uint32_t countAlpha_cigar2 = strtoul(cigar2.substr(0, digit_end_cigar2).c_str(),0,10);

	
	if (correspondTable[cigar[0]&0xf] ==cigar1[len_cigar1-1]) {
		cigar[0] += (countAlpha_cigar1 << 4);
		cigar1.erase(digit_start_cigar1+1,len_digit + 1);
	}

	if ( correspondTable[cigar[n_cigar-1]&0xf] == cigar2[digit_end_cigar2]) {
		cigar[n_cigar-1] += (countAlpha_cigar2 << 4) ;

		cigar2.erase( 0, digit_end_cigar2 + 1);
	}

	return 0;
}



int Aligner::connect(SvSam_Rec *rec1, SvSam_Rec *rec2, kseq_t *trunk)
{
	// cout<<"sdfasdfds"<<endl;
	//rec1->ref_end = rec2->ref_end;
	//rec1->read_end = rec2->read_end;
	uint8_t 	readqry[1024];
	uint8_t 	refqry[1024];
	string transitionPart = "";
	char 		tempCigar[50];
	int 		bscore = 0;
	char 		*readStartP;
	char 		*refStartP;

	int read_len = rec2->read_start - rec1->read_end;
	//cout<<rec2->read_start<<"\t"<<rec1->read_end<<endl;
	//cout<<rec2->ref_start<<"\t"<<rec1->ref_end<<endl;
	int ref_len = rec2->ref_start - rec1->ref_end;
	// fprintf(stderr,"%hu\t%hu\n",rec2->read_start,rec1->read_end);
	// fprintf(stderr,"%hu\t%hu\n",rec2->ref_start,rec1->ref_end);
	// fprintf(stderr,"%d\t%d\n",ref_len,read_len);
	uint32_t 	*cigar;
	int 		n_cigar = 0;

	const 	char 	correspondTable[] = "MIDNSHP=X";
	if (read_len == 0 || ref_len == 0) {
		if (read_len != 0) {
			//sprintf(tempCigar,"%dI",read_len);
			//transitionPart.append(tempCigar);
			//
			cigar = (uint32_t *)calloc(1,sizeof(uint32_t));
			cigar[0] = (read_len<<4) | 2;
			n_cigar = 1;
			bscore = opt->gapopen + (read_len - 1)*opt->gapextend;
		} else {
			if (0 != ref_len) {
				cigar = (uint32_t *)calloc(1,sizeof(uint32_t));
				cigar[0] = (ref_len<<4) | 3;
				n_cigar = 1;
				bscore = opt->gapopen + (ref_len - 1)*opt->gapextend;
			}
		}
	} else {

		if (rec1->flag)
			readStartP = trunk->seq.rs + rec1->read_end;
		else
			readStartP = trunk->seq.s + rec1->read_end;
		refStartP = genome + rec1->ref_start;
		transIntoDec(readqry,readStartP,read_len);
		transIntoDec(refqry,refStartP,ref_len);

		const uint8_t *readqry_ = readqry;
		const uint8_t *refqry_ = refqry;
		//prseq(readStartP,read_len,true);
		//prseq(refStartP,ref_len,true);
		
		int w = read_len > ref_len ? read_len : ref_len;
		//opt->gap open according to address method it may be slow, it will be changed soon
		bscore = ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,opt->gapopen,opt->gapextend,w,&n_cigar,&cigar);

		//fprintf(stderr,"%d %d %d %d %d %d\n",order[i-1],order[i],read_len,ref_len,score, n_cigar);
		//cout<<temp<<"\t"<<score<<'\t'<<"4"<<endl;
		
	}


	
	//merge those same alpha(s)

	//in case both ref_len and read_len are zero
 	if (n_cigar) {
 		mergeBoundaryCigar(rec1->cigar,rec2->cigar,cigar,n_cigar, correspondTable);

 		for (int z=0;z<n_cigar;++z) {
 		 	sprintf(tempCigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
 		// 	//sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
 		// 	//++startPosCigar;
 		 	transitionPart.append(tempCigar);
 		}		
 		free(cigar);	
 	}
	


	rec1->cigar += transitionPart + rec2->cigar;
	//rec1->tailCigar = rec2->tailCigar;
	//rec1->bodyScore += bscore + rec2->bodyScore;
	//rec1->tailScore = rec2->tailScore;
	rec1->score += bscore + rec2->score;
	rec1->ref_end = rec2->ref_end;
	rec1->read_end = rec2->read_end;


}
//input parameters char * read, char * rcRead for multiple threads

int Aligner::produceSAM(SvSam_Rec *_svsams , int countbulks,int *sam4bulk, kseq_t *trunk, uint32_t *len)
{
	int quality = 0;
	//bool besubstr[sam4bulk[countbulks]];
	int extend[sam4bulk[countbulks]];

	for (int i=0;i<sam4bulk[countbulks];++i) {
		extend[i] = 0;
		_svsams[i].MAQ = 0;
	}

	for (int i=0;i<countbulks;++i) {
		for (int j=sam4bulk[i];j<sam4bulk[i+1];++j) {
			if (~_svsams[j].MAQ) {
				for (int m=i+1;m<countbulks;++m) {
					int q;
					if ((q=qualified(_svsams+j,_svsams,sam4bulk[m],sam4bulk[m+1]))!= -1) {
						if (_svsams[j].flag) {
							connect(_svsams+q,_svsams+j,trunk);
							extend[q] = extend[j] - 1;
							_svsams[j].MAQ = -1;

						} else {
							connect(_svsams+j, _svsams+q, trunk);
							++extend[j];// something wrong?
							_svsams[q].MAQ = -1;
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
		if (~_svsams[i].MAQ) ++countOutput;
	}

	int secbestscore,secbestInd,secbestbulk;
	int bestscore,bestInd,bestbulk;
	bestscore = 0x80000001;

	for (int i=0;i<countbulks;++i) {
		for (int j=sam4bulk[i];j<sam4bulk[i+1];++j)
		if (~_svsams[j].MAQ) {
			if (bestscore <= _svsams[j].score) {
				secbestscore = bestscore;
				secbestInd = bestInd;
				secbestbulk = bestbulk;
				bestscore = _svsams[j].score;
				bestInd = j;
				bestbulk = i;

			} else {
				if (secbestscore <= _svsams[j].score) {
					secbestscore = _svsams[j].score;
					secbestInd = j;
					secbestbulk = i;
				}
			}
		}
	}
	//quality = 0;
//output first:
	//int rclip,lclip;
	//char *usedread;
	//char *usedqual;
	//char *probQual = "*";
	if (_svsams[bestInd].flag) {
		_svsams[bestInd].rclip = len[bestbulk+extend[bestInd]];
		_svsams[bestInd].lclip = trunk->seq.l - len[bestbulk+1];
		//usedread = trunk->seq;
	} else {
		_svsams[bestInd].rclip = len[bestbulk];
		_svsams[bestInd].lclip = trunk->seq.l - len[bestbulk+1+extend[bestInd]];
		//usedread = trunk->seq.s;
	}

	//if (trunk->qual.s == NULL)
	//	usedqual = probQual;
	//else
	//	usedqual = trunk->qual.s;

/*
	cout<<trunk->name.s<<"\t"<<_svsams[bestInd].flag<<"\t"<<ChrName[_svsams[bestInd].chrIndex]<<"\t"<<_svsams[bestInd].pos<<"\t"<<quality<<"\t";

	if (lclip)  	cout<<lclip<<"S";
	cout<<_svsams[bestInd].headCigar<<_svsams[bestInd].bodyCigar<<_svsams[bestInd].tailCigar;
	if (rclip) 	cout<<rclip<<"S";
	cout<<"\t"<<"*"<<"\t"<<"0"<<"\t"<<"0"<<"\t";

	cout<<usedqual<<"\t"<<"AS:i:"<<_svsams[bestInd].score<<endl;
*/
	if (countOutput!=1) {//
		for (int i=0;i<countbulks;++i) {
			for (int j=sam4bulk[i];j<sam4bulk[i+1];++j) {
				if (~_svsams[j].MAQ && j != bestInd) {
					//output
					if (_svsams[j].flag) {
						_svsams[j].rclip = len[i+extend[j]];
						_svsams[j].lclip = trunk->seq.l - len[i+1];
						//usedread = RCRead;
					} else {
						_svsams[j].lclip = len[i];
						_svsams[j].rclip = trunk->seq.l - len[i+1+extend[j]];
						//usedread = trunk->seq.s;
					}
					_svsams[j].flag += 256;

					// cout<<trunk->name.s<<"\t"<<_svsams[j].flag + 256<<"\t"<<ChrName[_svsams[j].chrIndex]<<"\t"<<_svsams[j].pos<<"\t"<<quality<<"\t";
					// if (lclip)	cout<<lclip<<"S";
					// cout<<_svsams[j].headCigar<<_svsams[j].bodyCigar<<_svsams[j].tailCigar;
					// if (rclip)	cout<<rclip<<"S";
					// cout<<"\t"<<"*"<<"\t"<<"0"<<"\t"<<"0"<<"\t";
					// cout<<usedread<<"\t";
					// cout<<usedqual<<"\t"<<"AS:i:"<<_svsams[j].score<<endl;
				}
			}
		}
	}
	//exchange pos 1 with pos best
	SvSam_Rec temp;
	temp = _svsams[bestInd];
	_svsams[bestInd] = _svsams[0];
	_svsams[0] = temp;
	return 0;
}

int Aligner::rhat_seq_read(kstream_t *_fp, kseq_t *_seqs, int n_needed)
{
	kseq_t *temp = _seqs;

	int i = 0;
	//temp[i].f = _fp;
	while(i <n_needed && (temp[i].f = _fp) && kseq_read(temp+i)>=0 ) {
		//int z = 0;
		//fprintf(stderr, "%s\t%d\n",temp[i].name.s, i);
		revComRead(temp[i].seq.s, temp[i].seq.rs, temp[i].seq.l);
		++i;
	}
	return i;

}

thread_aux *thread_initiate(int n_thread, RHashtable **rhashtab, RHashtable **rrhashtab, uint32_t *sed_rec, uint16_t *sed_hit_times,
	uint16_t *unused_bkt, opts *_opt, Aligner *aln)
{
	thread_aux *aux = new thread_aux[n_thread];

	for (int i= 0; i<n_thread; ++i) {
		aux[i].aln = aln;
		aux[i].opt = _opt;
		aux[i].com_var.sed_hit_times = sed_hit_times + i * LEN_BASES;
		aux[i].com_var.unused_bkt = unused_bkt + i * LEN_BASES;
		aux[i].com_var.sed_rec = sed_rec + i * LEN_BASES;
		aux[i].com_var.rhashtab = rhashtab[i];
		aux[i].com_var.rrhashtab = rrhashtab[i];
	}
	return aux;
}

int Aligner::applySV(kseq_t *trunk, RHashtable *rhashtab, RHashtable *rrhashtab, SvSam_Rec *_svsams, uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt)
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
	//uint32_t bulk_offset = groupNum * opt->len_limit;

	bkt_index[0] = 0;
	len[0] = 0;
	for (uint32_t i=0;i<pNumber -1;++i) len[i+1] = len[i] + LEN_BASES;
	len[pNumber] = len[pNumber-1] + leftLen;

	char *useread;
	bool isRC;
	int  movement ;

	for (uint32_t i = 0; i < pNumber - 1;++i) {
		useread = trunk->seq.s;
		isRC = false;
		movement = i*LEN_BASES;
		for (int j=0;j<=1;++j) {
			proCans(useread+movement,LEN_BASES, isRC,cansHeap, _sed_rec, _sed_hit_times, _unused_bkt);
			//useread = RCRead + bulk_offset;
			useread = trunk->seq.rs;
			movement = trunk->seq.l - len[i+1];
			isRC = true;
		}
		int samCounter = conductAlign(trunk,trunk->seq.s+i*LEN_BASES, trunk->seq.rs + movement, LEN_BASES, rhashtab,rrhashtab,cansHeap,_svsams+bkt_index[i]);
		//update svsams
		bkt_index[i+1] = bkt_index[i] + samCounter;
		for (int k=bkt_index[i];k<bkt_index[i+1];++k) {
			if (_svsams[k].flag) {
				 _svsams[k].read_start += movement; // here remember the read length can only be less than 65535 otherwise the structure has to be changed.
				 _svsams[k].read_end += movement;

			} else {
				_svsams[k].read_start += len[i];
				_svsams[k].read_end += len[i];
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
		proCans(useread+ movement + offset, LEN_BASES, isRC,cansHeap, _sed_rec, _sed_hit_times, _unused_bkt);
		useread = trunk->seq.rs ;
		movement = 0;
		isRC = true;
	}

	//len[pNumber] = len[pNumber - 1] + leftLen;
	bkt_index[pNumber] = bkt_index[pNumber-1] + conductAlign(trunk,trunk->seq.s+(pNumber-1)*LEN_BASES, trunk->seq.rs, leftLen,rhashtab,rrhashtab,cansHeap, _svsams + bkt_index[pNumber-1]);
	for (int k=bkt_index[pNumber-1];k<bkt_index[pNumber];++k) {
		if (_svsams[k].flag == 0) {
			_svsams[k].read_start += len[pNumber-1];
			_svsams[k].read_end += len[pNumber-1];
		}
	}
	if (bkt_index[pNumber] <= 0) return 0;
	//produceSAM(SvSam_Rec *_svsams , int countbulks,int *sam4bulk, kseq_t *trunk, uint32_t *len)
	produceSAM(_svsams, pNumber, bkt_index, trunk, len);
	return bkt_index[pNumber];

}

int  Aligner::applyNonSV(kseq_t *trunk,  RHashtable *rhashtab, RHashtable *rrhashtab, Sam_Rec *_sams,
	uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt)
{

	std::priority_queue<bkt2> 	cansHeap;

	uint32_t diff = trunk->seq.l > LEN_BASES?trunk->seq.l - LEN_BASES:0;

	uint32_t offset = diff>>1;

	uint32_t canReadLen = trunk->seq.l > LEN_BASES ? LEN_BASES:trunk->seq.l;

	char * useread = trunk->seq.s;
	bool isRC = false;
	//uint32_t bulk_offset = groupNum * opt->len_limit;
	//uint32_t extend = (uint32_t)(1.1*offset);//waiting to be tested;

	for (int j=0;j<=1;j++) {
		proCans(useread+offset, canReadLen, isRC, cansHeap, _sed_rec, _sed_hit_times, _unused_bkt);
		useread = trunk->seq.rs;
		isRC = true;
	}
	return  conductAlign(trunk,cansHeap, rhashtab,rrhashtab, _sams);
	//sign = trunk->seq.l < LEN_BASES ? 1:sign; //if read length is less than len bases even though it was not mapped it won't conduct sv operation;

}


static void *thread_worker(void *data)
{
	thread_aux *aux = (thread_aux *)data;
	int _read_seq;
	//for (int i=0;i<aux->n_seqs;++i) {
			//printf("%s\t%u\n", aux->seqs[i].name.s,aux->sam_details[i]);
	//}
	while (1) {
		pthread_rwlock_wrlock(&rwlock);
		_read_seq = read_seq++;
		//cout<<aux->tid<<"\t"<<_read_seq<<endl;
		pthread_rwlock_unlock(&rwlock);
		if (_read_seq < aux->n_seqs) {
			
			aux->sam_details[_read_seq] = 0;
			aux->sam_details[_read_seq] = aux->aln->applyNonSV(aux->seqs+_read_seq,aux->com_var.rhashtab, aux->com_var.rrhashtab, aux->sams + _read_seq * aux->opt->canN, aux->com_var.sed_rec, aux->com_var.sed_hit_times, aux->com_var.unused_bkt);
			if (!aux->sam_details[_read_seq] && !(aux->seqs[_read_seq].seq.l < LEN_BASES)) {
				if (!aux->svsams[_read_seq])
					aux->svsams[_read_seq] = new SvSam_Rec[aux->opt->canN *(aux->opt->len_limit/LEN_BASES + 1)];
					aux->sam_details[_read_seq] = aux->aln->applySV(aux->seqs+_read_seq,aux->com_var.rhashtab, aux->com_var.rrhashtab, aux->svsams[_read_seq], aux->com_var.sed_rec, aux->com_var.sed_hit_times, aux->com_var.unused_bkt);
				if (aux->sam_details[_read_seq])
					aux->sam_details[_read_seq] = aux->sam_details[_read_seq] << 1;
			} else
				aux->sam_details[_read_seq] = aux->sam_details[_read_seq] << 1 | 1;		
		} else break;	
	}
	//return 0;
	//return ;
}



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
	genome_e = genome + len_genome - 1;
	//output header
	cout<<"@HD\tVN:"<<PACKAGE_VERSION<<endl;
	for (int i=1;i<countChr;++i) {cout<<"@SQ\tSN:"<<ChrName[i]<<"\tLN:"<<Start_pos[i]-Start_pos[i-1]<<endl;}
	cout<<"@PG\tID:"<<PACKAGE_NAME<<"\tVN:"<<PACKAGE_VERSION<<"\tCL:";
	for (int i=0;i<opt->argc;++i) {cout<<opt->argv[i]<<" ";}
	cout<<endl;

	Hash hashh;
	hashtab = hashh.load_hashfile(opt->hashdir,len_genome,opt->len_sed);

	if ( NULL == hashtab) {
		fprintf(stderr,"Fail to load Hash Table, now exit");
		exit(1);
	}

	gzFile fp;
	//kseq_t *trunk;
	fp = gzopen(opt->readpath, "r");

	if (opt->thread <= 1) {

		kseq_t *seqs = kseq_init(fp);

		Sam_Rec *sams = new Sam_Rec[opt->canN];

		SvSam_Rec *svsams = new SvSam_Rec[opt->canN *(opt->len_limit/LEN_BASES + 1)];

		RHashtable *rhashtab = new RHashtable(opt->rh_seed_len,FORELEN,opt->len_limit);

		if ( NULL==rhashtab ) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		RHashtable *rrhashtab = new RHashtable(opt->rh_seed_len,FORELEN,opt->len_limit);

		if ( NULL==rrhashtab ) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		while (kseq_read(seqs)>=0) {

			revComRead(seqs->seq.s, seqs->seq.rs, seqs->seq.l);

			uint16_t sam_details = applyNonSV(seqs, rhashtab, rrhashtab, sams, sed_rec, sed_hit_times, unused_bkt);

			if (!sam_details && !(seqs->seq.l < LEN_BASES)) {

				sam_details = applySV(seqs, rhashtab, rrhashtab, svsams, sed_rec, sed_hit_times, unused_bkt);

				if (sam_details) sam_details = sam_details << 1;

			} else
			 	sam_details = sam_details << 1|1;
			//printf("%s\t%u\n",seqs->name.s,sam_details);
			outputSam(seqs, sams, &svsams, &sam_details, 1);
		}
		kseq_destroy(seqs);
	} else {
		kstream_t *_fp = ks_init(fp);

		int n_seqs;

		int n_needed = N_NEEDED;

		kseq_t 		*seqs = (kseq_t *)calloc(n_needed, sizeof(kseq_t));

		if (seqs == NULL)  { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		Sam_Rec 	*sams = new Sam_Rec[n_needed * opt->canN];

		if (sams == NULL)  { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		SvSam_Rec 	**svsams = new SvSam_Rec*[n_needed];

		if (svsams == NULL)  { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		for (int i=0;i<n_needed;++i) svsams[i] = NULL;

		uint16_t 	*sam_details = new uint16_t [n_needed];

		if (sam_details == NULL)  { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}

		RHashtable *rhashtab[opt->thread];

		for (int i=0;i<opt->thread;++i) {
			rhashtab[i] = new RHashtable(opt->rh_seed_len,FORELEN,opt->len_limit);
			if ( NULL==rhashtab[i] ) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}
		}

		RHashtable *rrhashtab[opt->thread];

		for (int i=0;i<opt->thread;++i) {
			rrhashtab[i] = new RHashtable(opt->rh_seed_len,FORELEN,opt->len_limit);
			if ( NULL== rrhashtab[i]) { fprintf(stderr, "Failed when applying for new space! now exit"); exit(1);}
		}


		thread_aux *aux = thread_initiate(opt->thread, rhashtab, rrhashtab, sed_rec, sed_hit_times, unused_bkt, opt, this);
		pthread_rwlock_init(&rwlock, NULL);
		while((n_seqs = rhat_seq_read(_fp, seqs, n_needed))>0) {
			read_seq = 0;
			pthread_t *tid;
			//pthread_attr_t attr;

			//pthread_attr_init(&attr);
			//pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			//data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
			tid = (pthread_t*)calloc(opt->thread, sizeof(pthread_t));
			//int average = n_seqs/opt->thread;
			//int lastOne = n_seqs - average * opt->thread + average;
			int j;
			for (j = 0; j < opt->thread; ++j) {
				aux[j].tid = j;
				aux[j].n_seqs = n_seqs;
				aux[j].seqs = seqs;
				aux[j].sams = sams;
				aux[j].svsams = svsams;
				aux[j].sam_details = sam_details;

				pthread_create(&tid[j], NULL, thread_worker, aux + j);
			}
			/*
			aux[j].tid = j;
			aux[j].n_seqs = lastOne;
			aux[j].seqs = seqs + j * average;
			aux[j].sams = sams + j * average * opt->canN;
			aux[j].svsams = svsams + j * average;
			aux[j].sam_details = sam_details + j * average;
	
			pthread_create(&tid[j], &attr, thread_worker, aux + j);
			*/
			for (j = 0; j < opt->thread; ++j) pthread_join(tid[j], 0);
			//free(data);
			free(tid);

			outputSam(seqs, sams, svsams, sam_details, n_seqs);
			//fprintf(stderr, "nothing wrong");

		}
		pthread_rwlock_destroy(&rwlock);
		free(seqs);
		delete []sams;
		delete []svsams;
		delete []sam_details;
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

	// delete rhashtab;
	// delete rrhashtab;

	gzclose(fp);


}
