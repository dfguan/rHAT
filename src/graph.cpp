/*
	kinda residule, have no Idea to solve it now.
*/
#include "graph.h"

#include <algorithm>
#include <iostream>
#include <cstdio>

using namespace std;
//#include "ksw.h"

#define WAITINGLEN 203
#define WAITINGLENLIMIT 437
#define STD_EXTRACT 400

#define ACCEPT_VERTEX 20
#define VERTEX_LIMIT 10000

//forelen should be less than 7 and bigger than 5

const uint8_t seq_nt4_tablet[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
//actually it could be neat, consider to remove it.
uint8_t transTable[] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 0,  4, 4, 4, 0,  4, 4, 4, 4,  4, 4, 1, 4,
	4, 4, 4, 4,  0, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 0,  4, 4, 4, 0,  4, 4, 4, 4,  4, 4, 1, 4,
	4, 4, 4, 4,  0, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


void prseq(char *str, int len, bool enter)
{
	for (int i=0;i<len;++i)  {
		cout<<str[i];
	}
	if (enter)
		cout<<endl;

}
RHashtable::RHashtable(uint32_t kmer,uint32_t forelen, uint32_t len_limit)
{
	limitOfp2leftSeq = (1<<(forelen << 1)) + 1;

	p2leftSeq = new uint16_t[limitOfp2leftSeq];

	left_seq = new uint16_t [len_limit];
	p2seqNum = new uint16_t [len_limit];
	seq_num = new uint16_t [len_limit];
	kmer_value = new uint32_t [len_limit];
	//seq_number = new uint16_t [len_limit];
	//order = new uint16_t [len_limit];
	len_sed = kmer;
	forelength = forelen;

}

RHashtable::~RHashtable()
{
	if (NULL != p2leftSeq) {
		delete [] p2leftSeq;
	}

	if (NULL != left_seq) {
		delete [] left_seq;
	}

	if (NULL != p2seqNum) {
		delete [] p2seqNum;
	}

	if (NULL != seq_num) {
		delete [] seq_num;
	}

	if (NULL != kmer_value) {
		delete [] kmer_value;
	}

}

int compare_kmer_value(const void *p, const void *q, void *t)
{
	uint16_t f = *(uint16_t *)p;
	uint16_t h = *(uint16_t *)q;
	uint32_t *s = (uint32_t *)t;
 	if (s[f] == s[h]) 	return f > h;
 	return s[f] - s[h];
}



void RHashtable::buildRHash(char *seq, uint32_t lenSeq)
{
	uint32_t pre = 0;
	uint32_t mask = 0xffffffff>>(32 - (len_sed<<1));



	//produce all kmer value;
	transfer(seq,&pre,len_sed);
	//prseq(seq,13,false);
	//cout<<"\t"<<0<<endl;
	kmer_value[0] = pre;

	//order[0] = 0;
	seq_num[0] = 0;
	// if you want to mask uncomment the next sentence
	for (uint32_t i=1; i<= lenSeq - len_sed; ++i) {
		pre = ((pre<<2)&mask)|trans[seq[i + len_sed - 1 ]];
		//prseq(i+seq, 13,false);
		//cout<<"\t"<<i<<endl;
		kmer_value[i] = pre;
		seq_num[i] = i;
		// int count_kmer

		//order[i] = i;
	}

	//sort them
	qsort_r(seq_num,lenSeq - len_sed + 1, sizeof(uint16_t), compare_kmer_value, kmer_value);

	// push them into rhash table
	uint32_t move = (len_sed - forelength) << 1;
	uint32_t temp = 0;
	uint32_t mask2 = (~temp) >> (32 - move);

	uint32_t pre_array_seq = kmer_value[seq_num[0]] >> move;
	uint16_t pre_ltseq_Value = (uint16_t)(kmer_value[seq_num[0]] & mask2);


	//put the first element into rhash table

	p2leftSeq[pre_array_seq] = 1;
	left_seq[0] = pre_ltseq_Value;
	p2seqNum[0] = 0;
	//seq_num[0] = order[0];

	uint32_t next_array_seq;
	uint16_t next_ltseq_value;
	uint32_t count_left_seq = 2;

	//initiate p2leftSeq
	for (uint32_t i=0;i<pre_array_seq;++i)
		p2leftSeq[i] = 0;


	for (uint32_t i=pre_array_seq+1;i<limitOfp2leftSeq;++i)
		p2leftSeq[i] = 0;

	for (uint32_t i = 1; i<= lenSeq - len_sed; ++i) {
		next_array_seq = kmer_value[seq_num[i]] >> move;
		next_ltseq_value = (uint16_t) (kmer_value[seq_num[i]]&mask2);

		if (next_array_seq == pre_array_seq) {
			if (next_ltseq_value != pre_ltseq_Value) {

				left_seq[count_left_seq - 1] = next_ltseq_value;
				p2seqNum[count_left_seq - 1] = i;

				++count_left_seq;
				pre_ltseq_Value = next_ltseq_value;
			}

		} else {

			p2leftSeq[pre_array_seq + 1] = count_left_seq;

			p2leftSeq[next_array_seq] = count_left_seq;
			left_seq[count_left_seq - 1] = next_ltseq_value;
			p2seqNum[count_left_seq - 1] = i;


			pre_array_seq = next_array_seq;
			pre_ltseq_Value = next_ltseq_value;
			++count_left_seq;
		}
		//seq_num[i] = order[i];
	}
	p2leftSeq[pre_array_seq + 1] = count_left_seq;
	p2seqNum[count_left_seq - 1] = lenSeq - len_sed + 1;

}

uint16_t binsearch(uint16_t sval,uint16_t low,uint16_t high,uint16_t *bkt)
{
  if (sval<bkt[low]||sval>bkt[high])
  	return 0xffff;
  while(low<=high)
  {
    uint32_t middle=low+((high-low)>>1);

    if(bkt[middle]<sval)
      low = middle + 1;
    else if(bkt[middle]>sval)
      high = middle - 1;
    else
      return middle;
  }
  return 0xffff;
}


uint16_t binsearchPos(uint16_t sval,uint16_t low,uint16_t high,uint16_t *bkt)
{
  if (sval<bkt[low])
  	return 0;
  else
  	if(sval>bkt[high])
  		return high;
  while(low<=high)
  {
    uint32_t middle=low+((high-low)>>1);

    if(bkt[middle]<sval)
      low = middle + 1;
    else if(bkt[middle]>sval)
      high = middle - 1;
    else
      return middle;
  }
  return low;
}
Graphic::Graphic(char *_ref_s, char *_ref_e) 
{
	ref_s = _ref_s;
	ref_t = _ref_e;
}
int Graphic::applyGraphic(RHashtable *rhashtab, char *ref, uint32_t lenRef, char *read, uint32_t lenRead,int *score, uint32_t waitingLen,
 uint32_t left_start,bool rc, uint32_t *startPos,  int countChr, Sam_Rec *sam, int countSam, int8_t *mat, int gapo, int gape)
{
	uint16_t seq_counter[lenRef];// = new uint16_t[lenRef];
	uint16_t p2startPos[lenRef]; //= new uint16_t [lenRef];
	for(uint32_t i=0;i<lenRef;++i) {
		seq_counter[i] = 0;
	}

	buildCounter(ref,lenRef,rhashtab,seq_counter,p2startPos);
	int flag;
	flag = createVertex(rhashtab->seq_num, ref, seq_counter,p2startPos, lenRef, 0,rhashtab->len_sed, read, lenRead,0, VERTEX_LIMIT);
	//cout<<flag<<endl;

	if (flag == -1) {
		if (lenRead < STD_EXTRACT) return -1;
		node.clear();
		livingSeed.clear();
		//cout<<"ccccc"<<endl;
		int sign;
		vertex head,tail;
		createVertex(rhashtab->seq_num, ref, seq_counter, p2startPos, lenRef - lenRead + STD_EXTRACT,0, rhashtab->len_sed,  read, STD_EXTRACT,0, -1);
		sign = findPos(lenRef - lenRead + STD_EXTRACT,STD_EXTRACT,waitingLen,true, &head);
		node.clear();
		livingSeed.clear();
		if (sign == -1)
			return -1;
		createVertex(rhashtab->seq_num, ref, seq_counter, p2startPos, lenRef - lenRead + STD_EXTRACT, lenRead - STD_EXTRACT, rhashtab->len_sed,  read, STD_EXTRACT,lenRead - STD_EXTRACT, -1);
		sign = findPos(lenRef - lenRead + STD_EXTRACT,STD_EXTRACT,waitingLen,false, &tail);
		tail.read_seq += lenRead - STD_EXTRACT;
		tail.ref_seq += lenRead - STD_EXTRACT;
		node.clear();
		livingSeed.clear();
		if (sign == -1)
			return -1;
		createLimVertex(rhashtab->seq_num,ref,lenRef,read,lenRead,seq_counter,p2startPos, rhashtab->len_sed, head,tail);
	}
	return dealGraph(lenRef, lenRead, read, ref, score, waitingLen, left_start, rc, startPos,  countChr, sam, countSam, mat, gapo, gape);
}

int Graphic::applyGraphic(RHashtable *rhashtab, char *ref, uint32_t lenRef, char *read, uint32_t lenRead,int *score, uint32_t waitingLen,
 uint32_t left_start,bool rc, uint32_t *startPos,  int countChr, SvSam_Rec *sam, int countSam, int8_t *mat, int gapo, int gape)
{
	uint16_t seq_counter[lenRef];// = new uint16_t[lenRef];
	uint16_t p2startPos[lenRef]; //= new uint16_t [lenRef];
	for(uint32_t i=0;i<lenRef;++i) {
		seq_counter[i] = 0;
	}

	buildCounter(ref,lenRef,rhashtab,seq_counter,p2startPos);
	int flag;
	flag = createVertex(rhashtab->seq_num, ref, seq_counter,p2startPos, lenRef, 0,rhashtab->len_sed, read, lenRead,0, VERTEX_LIMIT);
	//cout<<flag<<endl;

	if (flag == -1) {
		if (lenRead < STD_EXTRACT) return -1;
		node.clear();
		livingSeed.clear();
		//cout<<"ccccc"<<endl;
		int sign;
		vertex head,tail;
		createVertex(rhashtab->seq_num, ref, seq_counter, p2startPos, lenRef - lenRead + STD_EXTRACT,0, rhashtab->len_sed,  read, STD_EXTRACT,0, -1);
		sign = findPos(lenRef - lenRead + STD_EXTRACT,STD_EXTRACT,waitingLen,true, &head);
		node.clear();
		livingSeed.clear();
		if (sign == -1)
			return -1;
		createVertex(rhashtab->seq_num, ref, seq_counter, p2startPos, lenRef - lenRead + STD_EXTRACT, lenRead - STD_EXTRACT, rhashtab->len_sed,  read, STD_EXTRACT,lenRead - STD_EXTRACT, -1);
		sign = findPos(lenRef - lenRead + STD_EXTRACT,STD_EXTRACT,waitingLen,false, &tail);
		tail.read_seq += lenRead - STD_EXTRACT;
		tail.ref_seq += lenRead - STD_EXTRACT;
		node.clear();
		livingSeed.clear();
		if (sign == -1)
			return -1;
		createLimVertex(rhashtab->seq_num,ref,lenRef,read,lenRead,seq_counter,p2startPos, rhashtab->len_sed, head,tail);
	}
	return dealGraph(lenRef, lenRead, read, ref, score, waitingLen, left_start, rc, startPos,  countChr, sam, countSam, mat, gapo, gape);
}

void Graphic::buildCounter(char *seq, uint32_t len_seq, RHashtable *rhashtab,uint16_t *seq_counter, uint16_t *p2startPos)
{


	//prseq(seq,len_seq,true);
	uint32_t kmer = rhashtab->len_sed;
	uint32_t forelen = rhashtab->forelength;
	uint32_t mask = 0xffffffff>>(32 - (kmer<<1));


	uint32_t move = (kmer - forelen) << 1;
	uint32_t temp = 0;
	uint32_t mask2 = (~temp) >> (32 - move);




	uint16_t 	hitPos;
	//uint32_t seq_counter_count = 0;

	uint32_t 	pre;
	uint32_t 	array_seq;
	uint16_t 	leftValue;

	//prseq(seq,kmer,true);
	//first initiate

	uint8_t isIncludeN = 0 ;// the type of isIncludeN has to be changed if N

	for (uint32_t i=0;i<kmer;++i)  isIncludeN = (isIncludeN<<1)|transTable[seq[i]];

	transfer(seq,&pre,kmer);

	if (isIncludeN) seq_counter[0] = 0;
	else {
		array_seq = pre >> move;

		leftValue = (uint16_t)(pre & mask2);

		if (rhashtab->p2leftSeq[array_seq]&&rhashtab->p2leftSeq[array_seq+1]&&(rhashtab->p2leftSeq[array_seq]^rhashtab->p2leftSeq[array_seq+1])) {
			hitPos = binsearch(leftValue,rhashtab->p2leftSeq[array_seq]-1,rhashtab->p2leftSeq[array_seq+1] - 2,rhashtab->left_seq);
			if (hitPos^0xffff) {
				seq_counter[0] = rhashtab->p2seqNum[hitPos+1] - rhashtab->p2seqNum[hitPos];
				p2startPos[0] = rhashtab->p2seqNum[hitPos];
			}
		}
	}



	for (uint32_t i=1; i<= len_seq - kmer; ++i) {//check border
		//prseq(seq+i,kmer,true);
		pre = ((pre<<2)&mask)|trans[seq[i + kmer - 1 ]];
		isIncludeN = (isIncludeN<<1)|transTable[seq[i + kmer - 1]];
		if (isIncludeN) seq_counter[i] = 0;
		else {
			array_seq = pre >> move;
			leftValue = (uint16_t)(pre & mask2);
			//cout<<i<<"\t"<<rhashtab->p2leftSeq[array_seq]<<"\t"<<rhashtab->p2leftSeq[array_seq+1]<<endl;
			if (rhashtab->p2leftSeq[array_seq]&&rhashtab->p2leftSeq[array_seq+1]&&(rhashtab->p2leftSeq[array_seq]^rhashtab->p2leftSeq[array_seq+1])) {
				hitPos = binsearch(leftValue,rhashtab->p2leftSeq[array_seq]-1,rhashtab->p2leftSeq[array_seq+1] - 2,rhashtab->left_seq);
				//cout<<i<<"\t"<<hitPos<<endl;
				if (hitPos^0xffff) {
					seq_counter[i] = rhashtab->p2seqNum[hitPos+1] - rhashtab->p2seqNum[hitPos];
					p2startPos[i] = rhashtab->p2seqNum[hitPos];
				}
			}
		}
	}
}

void Graphic::createLimVertex(uint16_t *seq_n,char *ref, uint32_t lenRef, char *read, uint32_t lenRead,
	uint16_t *seq_counter,uint16_t *p2startPos, uint32_t kmer, vertex head, vertex tail)
{
	double indel = (double)(tail.read_seq - head.read_seq)/(double)(tail.ref_seq - head.ref_seq);

	vertex filled_element;
	ASeed  unrepeatElement;
	uint32_t counterOfLivingSeed;
	for (uint32_t i=0;i<=lenRef - kmer;++i) {
		//cout<<i<<"\t"<<seq_counter[i]<<endl;
		//for each living seed their life minus one point;
		//remove those point that has no life;
		//get how many point that own life
		// seq_counter[i] != Alive // use Alive to minus seq_counter[i]
		// binsearch (the seed find one ) not exist
		// extend it and add this one to living seed
		//else conitnue util test point is bigger than seq_counter wrong;
		//else you find

		counterOfLivingSeed = livingSeed.size();

		for (uint32_t w=0;w<counterOfLivingSeed;) {
			--livingSeed[w].left_time;
			if ( 0 == livingSeed[w].left_time) {
				livingSeed[w] = livingSeed[counterOfLivingSeed - 1];
				livingSeed.pop_back();
				--counterOfLivingSeed;
			}
			else {
				++livingSeed[w].read_seq;
				//++livingSeed[i].ref_seq;
				++w;
			}

		}
		//if (counterOfLivingSeed == 0) {
		//	cout<<"no seed for "<<i<<endl;

		//} else {
		//	for (uint32_t k=0;k<counterOfLivingSeed;++k) {
		//		cout<<livingSeed[k].read_seq<<"\t"<<livingSeed[k].left_time<<endl;
		//	}


		//}
		sort(livingSeed.begin(),livingSeed.end());

		//uint32_t test_counter = 0;

		if ( counterOfLivingSeed != seq_counter[i] ) {
			//whether it exist  in living seed:

			//if (seq_counter[i]<counterOfLivingSeed)
			// here I will use filter strategy.
			// find expected read pos:
			// binsearch its pos
			//reaturn all
			uint16_t expected_read_seq = head.read_seq + (uint16_t)(indel*(double)(i - head.ref_seq));

			uint32_t ind_beg = binsearchPos(expected_read_seq, 0, seq_counter[i]-1, seq_n + p2startPos[i]);


			uint32_t ind_end = seq_counter[i] <= ind_beg + ACCEPT_VERTEX - 1 ? seq_counter[i]:ind_beg + ACCEPT_VERTEX -1;

					 ind_beg = ind_beg < ACCEPT_VERTEX - 1 ? 0:ind_beg - ACCEPT_VERTEX + 1;

			uint32_t undiscovered = seq_counter[i] - counterOfLivingSeed;

			for(uint32_t j=ind_beg;j<ind_end;++j) {
				uint16_t seq = seq_n[p2startPos[i]+j];

				uint32_t t;
				for(t=0;t<counterOfLivingSeed;++t) {
					if (livingSeed[t].read_seq == seq)
						break;
				}
				if (t == counterOfLivingSeed) {

					uint32_t extend;

					for(extend=0; i+kmer+extend <lenRef && seq + kmer + extend < lenRead && ref[i+kmer+extend] == read[seq + kmer + extend];++extend);
					//prseq(ref+i,kmer+extend);
					//prseq(read+seq_n[p2startPos[i]+j],kmer+extend);
					filled_element.ref_seq = i;
					filled_element.read_seq = seq;
					filled_element.len = kmer + extend;
					node.push_back(filled_element);

					//unrepeatElement.ref_seq = i;
					unrepeatElement.read_seq = seq;
					unrepeatElement.left_time = extend + 1;
					livingSeed.push_back(unrepeatElement);
					//++vertex_counter;
					--undiscovered;
					if (undiscovered==0)
						break;
				}
			}
		}

	}
	return ;
}

int	Graphic::createVertex(uint16_t *seq_n, char *ref, uint16_t *seq_counter, uint16_t *p2startPos, uint32_t lenRef, uint32_t offset_ref,
	uint32_t kmer, char *read,uint32_t lenRead, uint32_t offset_read, uint32_t vertex_limit)
{
	vertex 		filled_element;
	ASeed  		unrepeatElement;
	uint32_t 	counterOfLivingSeed;
	//prseq(read,lenRead,true);
	uint32_t 	vertex_counter = 0;

	//read 		+= offset_read;
	ref 		+= offset_ref;
	seq_counter += offset_ref;
	p2startPos 	+= offset_ref;

	for (uint32_t i=0;i<=lenRef - kmer;++i) {
		//cout<<i<<"\t"<<seq_counter[i]<<endl;
		//for each living seed their life minus one point;
		//remove those point that has no life;
		//get how many point that own life
		// seq_counter[i] != Alive // use Alive to minus seq_counter[i]
		// binsearch (the seed find one ) not exist
		// extend it and add this one to living seed
		//else conitnue util test point is bigger than seq_counter wrong;
		//else you find
		if (vertex_counter > vertex_limit) {

			//cout<<vertex_counter<<endl;
			//node.clear();
			return -1;
		}
		counterOfLivingSeed = livingSeed.size();

		for (uint32_t w=0;w<counterOfLivingSeed;) {
			--livingSeed[w].left_time;
			if ( 0 == livingSeed[w].left_time) {
				livingSeed[w] = livingSeed[counterOfLivingSeed - 1];
				livingSeed.pop_back();
				--counterOfLivingSeed;
			}
			else {
				++livingSeed[w].read_seq;
				//++livingSeed[i].ref_seq;
				++w;
			}

		}
		//if (counterOfLivingSeed == 0) {
		//	cout<<"no seed for "<<i<<endl;

		//} else {
		//	for (uint32_t k=0;k<counterOfLivingSeed;++k) {
		//		cout<<livingSeed[k].read_seq<<"\t"<<livingSeed[k].left_time<<endl;
		//	}


		//}
		sort(livingSeed.begin(),livingSeed.end());

		//uint32_t test_counter = 0;
		if ( counterOfLivingSeed != seq_counter[i] ) {
			//whether it exist  in living seed:

			//if (seq_counter[i]<counterOfLivingSeed)

			uint32_t undiscovered = seq_counter[i] - counterOfLivingSeed;

			for(uint32_t j=0;j<seq_counter[i];++j) {
				uint16_t seq = seq_n[p2startPos[i]+j];

				if (seq < offset_read || seq > lenRead + offset_read - 1)
					continue;

				uint32_t t;
				for(t=0;t<counterOfLivingSeed;++t) {
					if (livingSeed[t].read_seq == seq)
						break;
				}
				if (t == counterOfLivingSeed) {

					uint32_t extend;
					//prseq(ref+i,kmer+1,true);
					//prseq(read+seq,kmer+1,true);
					for(extend=0; i+kmer+extend <lenRef && seq + kmer + extend < lenRead + offset_read && ref[i+kmer+extend] == read[seq + kmer + extend];++extend);
					//prseq(ref+i,kmer+extend);
					//prseq(read+seq_n[p2startPos[i]+j],kmer+extend);
					filled_element.ref_seq = i;
					filled_element.read_seq = seq - offset_read;
					filled_element.len = kmer + extend;
					node.push_back(filled_element);

					//unrepeatElement.ref_seq = i;
					unrepeatElement.read_seq = seq;
					unrepeatElement.left_time = extend + 1;
					livingSeed.push_back(unrepeatElement);
					++vertex_counter;
					--undiscovered;
					if (undiscovered==0)
						break;
				}
			}
		}

	}
	//exit(1);
	return 1;
}
void 	Graphic::revstr(uint8_t *revstring, char *string, int len)
{
	for (int i=0;i<len;++i) {
		revstring[i] = seq_nt4_tablet[string[len - 1 - i]];
	}

}
void 	Graphic::dealCigar(char *cigarbuf, char *headbuf, int headbuflen)
{
	int pre, next;
	next = headbuflen - 1;
	pre = next - 1;
	int posOfCigarBuf = 0;
	while(pre>=0) {
		while (pre>=0&&headbuf[pre]<='9') {
			--pre;
		}
		for (int z=pre+1;z<=next;z++) {
			cigarbuf[posOfCigarBuf++] = headbuf[z];
		}
		next = pre;
		pre = next - 1;
	}
	cigarbuf[posOfCigarBuf] = '\0';

}

int 	Graphic::transIntoDec(uint8_t *transtr,char *str, int length)
{
	for (int i=0;i<length;++i) {
		transtr[i] = seq_nt4_tablet[str[i]];
	}
	return 0;
}
int 	Graphic::CalEditDistancewithCigar(int *order, int order_len, char *read, uint32_t totalReadlen, char *ref, uint32_t totalReflen, uint32_t left_start,
		bool rc, uint32_t *chrStartP,  int countChr, SvSam_Rec *sams, int countSam, int8_t *mat, int gapo, int gape)
{
	//calculate 0 to order[order_len-1];
	//fprintf(stderr,"len%d",order_len);
	//int score 	= 	0 ;

	sams[countSam].cigar = "";
	//sams[countSam].tailCigar = "";
	//sams[countSam].headScore = 0;
	//sams[countSam].bodyScore = 0;
	sams[countSam].score = 0;


	char *readStartP;
	char *refStartP;
	//char *cigarStartP = cigarbuf;
	//int usedCigarsize;
	//int totalUsedCigarSize;
	//int refendpos = 0;
	int n_cigar = 0;
	uint32_t *cigar;
	int qlen = 0;
	int tlen = 0;
	//int gtle;
	//int max_off;
	//int gscore;
	//uint32_t endpos = node[order[1]].ref_seq + node[order[1]].len;
	uint32_t startpos;
	int w;
	uint32_t countM = 0;
	const 	char 	correspondTable[] = "MIDNSHP=X";
	const 	uint8_t *readqry_;
	const 	uint8_t *refqry_;

	int read_len 	= 	node[order[order_len-1]].read_seq;
	//figure out ref_len
	int ref_len 	= read_len;
	
	
	



	uint8_t  ind;
	uint32_t r_startP;
	uint32_t chrstartPos;
	char  	 trans_cigar[50];
	// 	 h0 = 0;
	if (read_len != 0) {


		//prseq(read,read_len,true);

		//char headbuf[(WAITINGLEN<<1)+1];
		revstr(readqry,read,read_len);
		//prseq(read,read_len,true);
		//transIntoDec(readqry,revreadqry,read_len);
		//refStartP = ref + node[order[order_len-1]].ref_seq - read_len;
		//prseq(refStartP,read_len,true);
		refStartP 		= 	ref + node[order[order_len-1]].ref_seq - read_len;
		if (refStartP < ref_s) { refStartP = ref_s; ref_len = left_start;}

		revstr(refqry,refStartP,ref_len);
		//cout<<node[order[order_len-1]].ref_seq;
		//prseq(refStartP,read_len,true);
		//transIntoDec(refqry,revrefqry,read_len);


		readqry_ = readqry;
		refqry_ = refqry;
		//h0 = read_len; 
		sams[countSam].score = ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
		//cout<<node[order[order_len-1]].ref_seq<<'\t'<<tlen<<"\t";
		//cout<<score<<endl;

	}
	startpos = node[order[order_len-1]].ref_seq - tlen;
	r_startP = left_start + startpos;
	for (ind=1;ind<countChr;++ind)
		if (chrStartP[ind]>r_startP)
			break;
	chrstartPos = r_startP - chrStartP[ind-1];
	// if exists;

	for (uint8_t i=0; i<countSam; ++i) {
		if (ind == sams[i].chrIndex && chrstartPos == sams[i].pos)
			return 0;
	}
	sams[countSam].chrIndex = ind;
	sams[countSam].pos = chrstartPos;
	sams[countSam].ref_start = chrstartPos;
	//sams[countSam].read_start = read_len;
	//cout<<seqN<<"\t";
	if (rc)
		sams[countSam].flag = 16;
	else
		sams[countSam].flag = 0;

		//cout<<16<<"\t"; else cout<<0<<"\t";
	//cout<<chrName[ind]<<"\t"<<chrstartPos<<"\t"<<0<<"\t";

	//for(int z=n_cigar-1;z>=0;--z)
	//	cout<<cigar[z];
	//cout<<endl;
	int startPosCigar = 0;
	sams[countSam].read_start = read_len - qlen;
	

	if (n_cigar) {
		for (int z=n_cigar-1;z>0;--z) {
			//cout<<(cigar[z]>>4)<<correspondTable[cigar[z]&0xf];
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
			//++startPosCigar;
		}

		if (correspondTable[cigar[0]&0xf] == 'M')
			countM = cigar[0] >> 4;
		else {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
		}

		free(cigar);
	}
	//if (n_cigar!=0) free(cigar);
	//deal with middle part

	for (int i= order_len-1;i>1;--i) {
		// deal with the same part
		//cigarStartP[0] = node[order[i]].len;
		//cigarStartP[1] = 'M'
		//stringLen = sprintf(cigarStartP,"%d",node[order[i]].len);
		//cigarStartP += stringLen;
		//cigarStartP[0] = 'M';

		//cigarbuflen -= (stringLen + 1);
		//cigarStartP += 1;
		//cout<<node[order[i]].len<<'M';
		
		countM += node[order[i]].len;

		//sams[countSam]._cigar[startPosCigar] = 'M';
		//++startPosCigar;
		//deal with different one
		sams[countSam].score += node[order[i]].len;
		//cout<<node[order[i]].len<<"\t"<<score<<'\t'<<"1"<<endl;
		readStartP = read + node[order[i]].read_seq + node[order[i]].len;
		refStartP = ref + node[order[i]].ref_seq + node[order[i]].len;
		read_len = node[order[i-1]].read_seq - (node[order[i]].read_seq + node[order[i]].len);
		ref_len = node[order[i-1]].ref_seq - (node[order[i]].ref_seq + node[order[i]].len);
		if ( 0 == read_len || 0 == ref_len) {
			//stick 'M' first
			if (read_len || ref_len) {
				startPosCigar += sprintf(trans_cigar, "%uM",countM);
				sams[countSam].cigar.append(trans_cigar);
				countM = 0;
				if (0 != read_len) {

					//stringLen = sprintf(cigarStartP,"%d",read_len);
					//cigarStartP += stringLen;
					//cigarStartP[0] = 'D';
					//cigarbuflen -= (stringLen + 1);
					//cigarStartP += 1;

					//score += read_len;
					//cout<<read_len<<'D';
					startPosCigar += sprintf(trans_cigar, "%uI",read_len);
					//sams[countSam]._cigar[startPosCigar] = 'D';
					//++startPosCigar;
					sams[countSam].cigar.append(trans_cigar);
					sams[countSam].score -= gapo + (read_len - 1)*gape;
					//cout<<"\t"<<score<<'\t'<<"2"<<endl;

				} else if (0 != ref_len) {
					startPosCigar += sprintf(trans_cigar, "%uD",ref_len);
					//sams[countSam]._cigar[startPosCigar] = 'I';
					//++startPosCigar;
					sams[countSam].cigar.append(trans_cigar);
					sams[countSam].score -= gapo + (ref_len - 1)*gape;
				}
			
			}

		} else {
			transIntoDec(readqry,readStartP,read_len);
			transIntoDec(refqry,refStartP,ref_len);

			readqry_ = readqry;
			refqry_ = refqry;
			//prseq(readStartP,read_len,true);
			//prseq(refStartP,ref_len,true);

			w = read_len > ref_len ? read_len : ref_len;

			sams[countSam].score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w,&n_cigar,&cigar);

			//fprintf(stderr,"%d %d %d %d %d %d\n",order[i-1],order[i],read_len,ref_len,score, n_cigar);
			//first cigar[0] is M
			if (n_cigar-1) {
				if (correspondTable[cigar[0]&0xf] == 'M') {
					countM += (cigar[0] >> 4);
					startPosCigar += sprintf(trans_cigar,"%uM",countM);
					sams[countSam].cigar.append(trans_cigar);
				} else {
					startPosCigar += sprintf(trans_cigar,"%uM",countM);
					sams[countSam].cigar.append(trans_cigar);
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
					sams[countSam].cigar.append(trans_cigar);
				}
				countM = 0;
				for (int z=1;z<n_cigar-1;++z) {
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
					//sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
					//++startPosCigar;
					sams[countSam].cigar.append(trans_cigar);
				}

				if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
					countM = cigar[n_cigar-1] >> 4;
					//startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
					//sams[countSam].cigar.append(trans_cigar);
				} else {
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
					sams[countSam].cigar.append(trans_cigar);
				}

			} else {
				if (correspondTable[cigar[0]&0xf] == 'M') 
					countM += (cigar[0] >> 4);
				else {
					startPosCigar += sprintf(trans_cigar,"%uM",countM);
					sams[countSam].cigar.append(trans_cigar);
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
					sams[countSam].cigar.append(trans_cigar);
					countM = 0;

				}
			}
			free(cigar);
			//cigarbuflen -= usedCigarsize;
			//cigarStartP += usedCigarsize;
		}

	}

	//stringLen = sprintf(cigarStartP,"%d",node[order[1]].len);
	//cigarStartP += stringLen;
	//cigarStartP[0] = 'M';//

	//cigarbuflen -= (stringLen + 1);
	//cigarStartP += 1;
	//cout<<node[order[1]].len<<'M';
	//startPosCigar += sprintf(trans_cigar, "%uM",node[order[1]].len);
	//sams[countSam]._cigar[startPosCigar] = 'M';
	//++startPosCigar;
	//sams[countSam].cigar.append(trans_cigar);
	countM += node[order[1]].len;

	sams[countSam].score += node[order[1]].len;
	//cout<<"\t"<<score<<'\t'<<"5"<<endl;
	readStartP = read + node[order[1]].read_seq + node[order[1]].len;
	refStartP = ref + node[order[1]].ref_seq + node[order[1]].len;

	read_len = totalReadlen - (node[order[1]].read_seq + node[order[1]].len);
	//may be discussed later
	qlen = 0;
	tlen = 0;
	n_cigar = 0;
	if (0 != read_len) { // if without else may be it will display previous cigar {
		transIntoDec(readqry,readStartP,read_len);

		ref_len = refStartP + read_len - 1 > ref_t? ref_t - refStartP + 1: read_len;

		transIntoDec(refqry,refStartP,ref_len);

		readqry_ = readqry;
		refqry_ = refqry;

		sams[countSam].score = ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
		//score += ksw_global(read_len,readStartP,ref_len,refStartP,5,mat,GAPOPEN,GAPEXTENDED,read_len,&n_cigar,&cigar);
		//fprintf(stderr,"%d\n",score);
		
		//endpos = node[order[1]].ref_seq + node[order[1]].len + tlen;
		//cout<<"\t"<<score<<'\t'<<"6"<<endl;
	} 
	//cout<<endl<<score<<endl;
	if (n_cigar) {
		if (correspondTable[cigar[0]&0xf] == 'M') {
			countM += (cigar[0] >> 4);
			startPosCigar += sprintf(trans_cigar,"%uM",countM);
			sams[countSam].cigar.append(trans_cigar);
		} else {
			startPosCigar += sprintf(trans_cigar,"%uM",countM);
			sams[countSam].cigar.append(trans_cigar);
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
		}
		for (int z=1;z<n_cigar;++z) {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
		}
		free(cigar);
	} else {
			startPosCigar += sprintf(trans_cigar,"%uM",countM);
			sams[countSam].cigar.append(trans_cigar);
	}
	//sams[countSam]._cigar[startPosCigar] = '\0';
	sams[countSam].ref_end = chrstartPos + node[order[1]].ref_seq + node[order[1]].len + tlen;
	sams[countSam].read_end = node[order[1]].read_seq + node[order[1]].len + qlen;
	//sams[countSam].score = sams[countSam].headScore + sams[countSam].bodyScore + sams[countSam].tailScore;
	//cout<<"\t"<<"*"<<"\t"<<0<<"\t"<<0<<"\t"<<"*"<<"\t"<<seqQual<<endl;
	//cout<<sams[countSam]._cigar<<endl;
	return 1;
}

int 	Graphic::CalEditDistancewithCigar(int *order, int order_len, char *read, uint32_t totalReadlen, char *ref, uint32_t totalReflen, uint32_t left_start,
		bool rc, uint32_t *chrStartP,  int countChr, Sam_Rec *sams, int countSam, int8_t *mat, int gapo, int gape)
{
	//calculate 0 to order[order_len-1];
	//fprintf(stderr,"len%d",order_len);
	//int score 	= 	0 ;

	sams[countSam].cigar = "";
	//sams[countSam].tailCigar = "";
	//sams[countSam].headScore = 0;
	//sams[countSam].bodyScore = 0;
	sams[countSam].score = 0;


	char *readStartP;
	char *refStartP;
	//char *cigarStartP = cigarbuf;
	//int usedCigarsize;
	//int totalUsedCigarSize;
	//int refendpos = 0;
	int n_cigar = 0;
	uint32_t *cigar;
	int qlen = 0;
	int tlen = 0;
	//int gtle;
	//int max_off;
	//int gscore;
	//uint32_t endpos = node[order[1]].ref_seq + node[order[1]].len;
	uint32_t startpos;
	int w;
	const 	char 	correspondTable[] = "MIDNSHP=X";
	const 	uint8_t *readqry_;
	const 	uint8_t *refqry_;
    uint32_t  countM = 0;
	int read_len 	= 	node[order[order_len-1]].read_seq;
	int ref_len		= 	read_len;
	


	uint8_t  ind;
	uint32_t r_startP;
	uint32_t chrstartPos;
	char  	 trans_cigar[50];
	int h0 = 0;
	
	if (read_len != 0 ) {

		//prseq(read,read_len,true);

		//char headbuf[(WAITINGLEN<<1)+1];
		revstr(readqry,read,read_len);
		//prseq(read,read_len,true);
		//transIntoDec(readqry,revreadqry,read_len);
		refStartP = ref + node[order[order_len-1]].ref_seq - read_len;
		if (refStartP < ref_s) { refStartP = ref_s; ref_len = left_start;}
		//prseq(refStartP,read_len,true);
		revstr(refqry,refStartP,ref_len);
		//cout<<node[order[order_len-1]].ref_seq;
		//prseq(refStartP,read_len,true);
		//transIntoDec(refqry,revrefqry,read_len);


		readqry_ = readqry;
		refqry_ = refqry;
		
		sams[countSam].score = ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len, &qlen, &tlen, &cigar, &n_cigar)  - read_len;
		//cout<<node[order[order_len-1]].ref_seq<<'\t'<<tlen<<"\t";
		//cout<<score<<endl;

	}
	startpos = node[order[order_len-1]].ref_seq - tlen;
	r_startP = left_start + startpos;
	for (ind=1;ind<countChr;++ind)
		if (chrStartP[ind]>r_startP)
			break;
	chrstartPos = r_startP - chrStartP[ind-1];
	// if exists;

	for (uint8_t i=0; i<countSam; ++i) {
		if (ind == sams[i].chrIndex && chrstartPos == sams[i].pos)
			return 0;
	}
	sams[countSam].chrIndex = ind;
	sams[countSam].pos = chrstartPos;
	//sams[countSam].ref_start = chrstartPos + tlen;
	//sams[countSam].read_start = read_len;
	//cout<<seqN<<"\t";
	if (rc)
		sams[countSam].flag = 16;
	else
		sams[countSam].flag = 0;

		//cout<<16<<"\t"; else cout<<0<<"\t";
	//cout<<lplihrName[ind]<<"\t"<<chrstartPos<<"\t"<<0<<"\t";

	//for(int z=n_cigar-1;z>=0;--z)
	//	cout<<cigar[z];
	//cout<<endl;
	int startPosCigar = 0;

	if (qlen != read_len ) {
		startPosCigar += sprintf(trans_cigar, "%uS", read_len - qlen);
		sams[countSam].cigar.append(trans_cigar);
	}// proves that softclipings do exist

	if (n_cigar) {
		for (int z=n_cigar-1;z>0;--z) {
			//cout<<(cigar[z]>>4)<<correspondTable[cigar[z]&0xf];
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
			//++startPosCigar;
		}

		if (correspondTable[cigar[0]&0xf] == 'M')
			countM = cigar[0] >> 4;
		else {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
		}

		free(cigar);
	}
	//if (n_cigar!=0) free(cigar);
	//deal with middle part

	for (int i= order_len-1;i>1;--i) {
		// deal with the same part
		//cigarStartP[0] = node[order[i]].len;
		//cigarStartP[1] = 'M'
		//stringLen = sprintf(cigarStartP,"%d",node[order[i]].len);
		//cigarStartP += stringLen;
		//cigarStartP[0] = 'M';

		//cigarbuflen -= (stringLen + 1);
		//cigarStartP += 1;
		//cout<<node[order[i]].len<<'M';
		
		countM += node[order[i]].len;

		//sams[countSam]._cigar[startPosCigar] = 'M';
		//++startPosCigar;
		//deal with different one
		sams[countSam].score += node[order[i]].len;
		//cout<<node[order[i]].len<<"\t"<<score<<'\t'<<"1"<<endl;
		readStartP = read + node[order[i]].read_seq + node[order[i]].len;
		refStartP = ref + node[order[i]].ref_seq + node[order[i]].len;
		read_len = node[order[i-1]].read_seq - (node[order[i]].read_seq + node[order[i]].len);
		ref_len = node[order[i-1]].ref_seq - (node[order[i]].ref_seq + node[order[i]].len);
		if ( 0 == read_len || 0 == ref_len) {
			//stick 'M' first
			//can't both be zero
			if (read_len || ref_len){
				
				startPosCigar += sprintf(trans_cigar, "%uM",countM);
				sams[countSam].cigar.append(trans_cigar);
				countM = 0;
				if (0 != read_len) {

					//stringLen = sprintf(cigarStartP,"%d",read_len);
					//cigarStartP += stringLen;
					//cigarStartP[0] = 'D';
					//cigarbuflen -= (stringLen + 1);
					//cigarStartP += 1;

					//score += read_len;
					//cout<<read_len<<'D';
					startPosCigar += sprintf(trans_cigar, "%uI",read_len);
					//sams[countSam]._cigar[startPosCigar] = 'D';
					//++startPosCigar;
					sams[countSam].cigar.append(trans_cigar);
					sams[countSam].score -= gapo + (read_len - 1)*gape;
					//cout<<"\t"<<score<<'\t'<<"2"<<endl;

				} else if (0 != ref_len) {
					startPosCigar += sprintf(trans_cigar, "%uD",ref_len);
					//sams[countSam]._cigar[startPosCigar] = 'I';
					//++startPosCigar;
					sams[countSam].cigar.append(trans_cigar);
					sams[countSam].score -= gapo + (ref_len - 1)*gape;
				}	
			
			}

		} else {
			transIntoDec(readqry,readStartP,read_len);
			transIntoDec(refqry,refStartP,ref_len);

			readqry_ = readqry;
			refqry_ = refqry;
			//prseq(readStartP,read_len,true);
			//prseq(refStartP,ref_len,true);

			w = read_len > ref_len ? read_len : ref_len;

			sams[countSam].score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w,&n_cigar,&cigar);

			//fprintf(stderr,"%d %d %d %d %d %d\n",order[i-1],order[i],read_len,ref_len,score, n_cigar);
			//first cigar[0] is M
			if (n_cigar-1) {
				if (correspondTable[cigar[0]&0xf] == 'M') {
					countM += (cigar[0] >> 4);
					startPosCigar += sprintf(trans_cigar,"%uM",countM);
					sams[countSam].cigar.append(trans_cigar);
				} else {
					startPosCigar += sprintf(trans_cigar,"%uM",countM);
					sams[countSam].cigar.append(trans_cigar);
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
					sams[countSam].cigar.append(trans_cigar);
				}
				countM = 0;
				for (int z=1;z<n_cigar-1;++z) {
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
					//sams[countSam]._cigar[startPosCigar] = correspondTable[cigar[z]&0xf];
					//++startPosCigar;
					sams[countSam].cigar.append(trans_cigar);
				}

				if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
					countM = cigar[n_cigar-1] >> 4;
					//startPosCigar += sprintf(trans_cigar,"%u%c",countM,'M');
					//sams[countSam].cigar.append(trans_cigar);
				} else {
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
					sams[countSam].cigar.append(trans_cigar);
				}

			} else {
				if (correspondTable[cigar[0]&0xf] == 'M') 
					countM += (cigar[0] >> 4);
				else {
					startPosCigar += sprintf(trans_cigar,"%uM",countM);
					sams[countSam].cigar.append(trans_cigar);
					startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
					sams[countSam].cigar.append(trans_cigar);
					countM = 0;

				}
			}
			free(cigar);
			//cigarbuflen -= usedCigarsize;
			//cigarStartP += usedCigarsize;
		}

	}

	//stringLen = sprintf(cigarStartP,"%d",node[order[1]].len);
	//cigarStartP += stringLen;
	//cigarStartP[0] = 'M';//

	//cigarbuflen -= (stringLen + 1);
	//cigarStartP += 1;
	//cout<<node[order[1]].len<<'M';
	//startPosCigar += sprintf(trans_cigar, "%uM",node[order[1]].len);
	//sams[countSam]._cigar[startPosCigar] = 'M';
	//++startPosCigar;
	//sams[countSam].cigar.append(trans_cigar);
	countM += node[order[1]].len;

	sams[countSam].score += node[order[1]].len;
	//cout<<"\t"<<score<<'\t'<<"5"<<endl;
	readStartP = read + node[order[1]].read_seq + node[order[1]].len;
	refStartP = ref + node[order[1]].ref_seq + node[order[1]].len;

	read_len = totalReadlen - (node[order[1]].read_seq + node[order[1]].len);
	//may be discussed later
	qlen = 0;
	tlen = 0;
	n_cigar = 0;
	if (0 != read_len) { // if without else may be it will display previous cigar {
		transIntoDec(readqry,readStartP,read_len);

		ref_len = refStartP + read_len - 1 > ref_t? ref_t - refStartP + 1: read_len;

		transIntoDec(refqry,refStartP,ref_len);

		readqry_ = readqry;
		refqry_ = refqry;

		sams[countSam].score = ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
		//score += ksw_global(read_len,readStartP,ref_len,refStartP,5,mat,GAPOPEN,GAPEXTENDED,read_len,&n_cigar,&cigar);
		//fprintf(stderr,"%d\n",score);
		
		//endpos = node[order[1]].ref_seq + node[order[1]].len + tlen;
		//cout<<"\t"<<score<<'\t'<<"6"<<endl;
	} 
	//cout<<endl<<score<<endl;
	if (n_cigar) {
		if (correspondTable[cigar[0]&0xf] == 'M') {
			countM += (cigar[0] >> 4);
			startPosCigar += sprintf(trans_cigar,"%uM",countM);
			sams[countSam].cigar.append(trans_cigar);
		} else {
			startPosCigar += sprintf(trans_cigar,"%uM",countM);
			sams[countSam].cigar.append(trans_cigar);
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
		}
		for (int z=1;z<n_cigar;++z) {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			sams[countSam].cigar.append(trans_cigar);
		}
		free(cigar);
	} else {
			startPosCigar += sprintf(trans_cigar,"%uM",countM);
			sams[countSam].cigar.append(trans_cigar);
	}

	if (qlen != read_len) {
		startPosCigar += sprintf(trans_cigar, "%uS", read_len - qlen);
		sams[countSam].cigar.append(trans_cigar);
	}
	//cout<<endl<<score<<endl;
	
	//sams[countSam]._cigar[startPosCigar] = '\0';
	//sams[countSam].ref_end = chrstartPos + node[order[1]].ref_seq + node[order[1]].len - startpos;
	//sams[countSam].read_end = node[order[1]].read_seq + node[order[1]].len;
	//sams[countSam].score = sams[countSam].headScore + sams[countSam].bodyScore + sams[countSam].tailScore;
	//cout<<"\t"<<"*"<<"\t"<<0<<"\t"<<0<<"\t"<<"*"<<"\t"<<seqQual<<endl;
	//cout<<sams[countSam]._cigar<<endl;
	return 1;
}

int 	Graphic::findPos(uint32_t lenRef, uint32_t lenRead, uint32_t waitingLen,bool type, vertex *vnode)
{
	vertex start,end;
	start.read_seq = 0;
	start.ref_seq = 0;
	start.len = 0;

	end.read_seq = lenRead - 1;
	end.ref_seq = lenRef - 1;
	end.len = 0;

	int counterP = 0;

	int nodeSize = node.size();

	sort(node.begin(),node.end());

	node.insert(node.begin(),start);//maybe changed
	node.push_back(end);


	//cout<<nodeSize<<endl;
	//prseq(read,lenRead,true);
	//prseq(ref,lenRef,true);

	// cout<<"After sorting:"<<endl;
	// for (int i=0;i<nodeSize;++i) {
	// 	cout<<node[i].read_seq<<"\t"<<node[i].ref_seq<<"\t"<<node[i].len<<endl;
	// }
	// cout<<endl;

	int   backtrace[nodeSize];
	short coverage[nodeSize];


	for(int i=0;i<nodeSize;++i) {
		//order[i] = i;
		backtrace[i] = -1;
		coverage[i] =  0;
	}

	//qsort_r(order,nodeSize,sizeof(uint32_t),compare_riAti,&node[0]);


	//node[i].read_seq  < WAITINGLEN + node[j].read_seq + node[j].len;
	//node[j].ref_seq + node[j].len < node[i].ref_seq;

	//connect to head
	for (int i=0;i<nodeSize;++i) {
		if (node[i].read_seq <= waitingLen) {
			backtrace[i] = 0;
			++counterP;
		}
		else
			break;
	}
/*
	if (0 == counterP) {
		for (int i=0;i<nodeSize;++i) {
			if (node[i].read_seq <= WAITINGLENLIMIT) {
				backtrace[i] = 0;
				++counterP;
			}
			else
				break;
		}
	}
*/
	if ( 0 == counterP ) return -1;
	//connect middle part
	for (int i=1; i< nodeSize - 1; i++) {
		//counterP = 0;
		for (int j= i-1;j>=1;--j) {
			if (node[i].read_seq >= node[j].read_seq + node[j].len) {
				if (node[i].read_seq <= waitingLen + node[j].read_seq + node[j].len) {
					if (node[j].ref_seq + node[j].len <= node[i].ref_seq && node[i].ref_seq <= node[j].ref_seq + node[j].len + waitingLen ) {
					//calculate coverage;
						//++counterP;
						if (coverage[i] < coverage[j] + node[j].len) {
							coverage[i] = coverage[j] + node[j].len;
							backtrace[i] = j;
							//cout<<backtrace[i]<<"->"<<i<<"\t"<<coverage[i]<<endl;
						}
					}

				} else 	break;
			}
		}
/*
		if (0==counterP) {
			for (int j= i-1;j>=0;--j) {
				if (node[i].read_seq >= node[j].read_seq + node[j].len) {
					if (node[i].read_seq <= WAITINGLENLIMIT + node[j].read_seq + node[j].len) {
						if (node[j].ref_seq + node[j].len <= node[i].ref_seq) {
							//calculate coverage;
							++counterP;
							if (coverage[i] < coverage[j] + node[j].len) {
								coverage[i] = coverage[j] + node[j].len;
								backtrace[i] = j;
								//cout<<backtrace[i]<<"->"<<i<<"\t"<<coverage[i]<<endl;
							}
						}

					} else 	break;
				}
			}
		}
*/
		//if (0 == counterP) return -1;

	}


	//connnect to tail
	counterP = 0;
	for (int i=nodeSize -2;i>=1;--i) {
		if (node[i].read_seq + waitingLen + node[i].len > node[nodeSize -1].read_seq ) {
			++counterP;
			if (coverage[nodeSize-1] < coverage[i] + node[i].len) {
				coverage[nodeSize-1] = coverage[i] + node[i].len;
				backtrace[nodeSize-1] = i;

			}
		} else break;
	}
/*
	if (0 == counterP) {
		for (int i=nodeSize -2;i>=0;--i) {
			if (node[i].read_seq + WAITINGLENLIMIT + node[i].len > node[nodeSize -1].read_seq) {
				++counterP;
				if (coverage[nodeSize-1] < coverage[i] + node[i].len) {
					coverage[nodeSize-1] = coverage[i] + node[i].len;
					backtrace[nodeSize-1] = i;

				}
			} else break;
		}
	}
*/

	if (0 == counterP) return -1;

	int 			rightOrder[nodeSize];
	int   			counterofRightOrder;

	//cout<<"backtrace: "<<nodeSize-1<<"->";
	rightOrder[0] = nodeSize - 1;
	counterofRightOrder = 1;
	for(int k= backtrace[nodeSize-1];k!=0;k=backtrace[k]) {
		if (k<0) {
			return -1;
			//break;
		} else
			rightOrder[counterofRightOrder++] = k;
			//cout<<k<<"->";
	}

	//for (int k=counterofRightOrder -1;k>=0;--k) {
	//	cout<<rightOrder[k]<<"->";
	//}
	//cout<<endl;
	//for (int i=counterofRightOrder-1;i>=1;--i)
	//	cout<<node[rightOrder[i]].read_seq<<"\t"<<node[rightOrder[i]].ref_seq<<endl;
	//fprintf(stderr,"%d",counterofRightOrder);
	if (counterofRightOrder <= 1)
		return -1;
	if (type) {
		vnode->read_seq =  node[rightOrder[counterofRightOrder-1]].read_seq;
		vnode->ref_seq = node[rightOrder[counterofRightOrder-1]].ref_seq;
		vnode->len = node[rightOrder[counterofRightOrder-1]].len;
	} else {
		vnode->read_seq =  node[rightOrder[1]].read_seq;
		vnode->ref_seq = node[rightOrder[1]].ref_seq;
		vnode->len = node[rightOrder[1]].len;
	}
	// ?????? this may be change after I amend this program ???? so just comment it here to alter myself.
	//after this action has been token, node should be cleared since it may be used by the other intention like figure out tail or head.
	node.clear();
	return 1;
}


int 	Graphic::dealGraph(uint32_t lenRef, uint32_t lenRead, char *read, char *ref, int *score, uint32_t waitingLen, uint32_t left_start,
		bool rc, uint32_t *startPos, int countChr, Sam_Rec *sam, int countSam, int8_t *matrix, int gapo, int gape)
{
	vertex start,end;
	start.read_seq = 0;
	start.ref_seq = 0;
	start.len = 0;

	end.read_seq = lenRead - 1;
	end.ref_seq = lenRef - 1;
	end.len = 0;

	int counterP = 0;


	sort(node.begin(),node.end());

	node.insert(node.begin(),start);//maybe changed
	node.push_back(end);

	int nodeSize = node.size();
	//cout<<nodeSize<<endl;
	//prseq(read,lenRead,true);
	//prseq(ref,lenRef,true);

	/*cout<<"After sorting:"<<endl;
	for (int i=0;i<nodeSize;++i) {
		cout<<node[i].read_seq<<"\t"<<node[i].ref_seq<<"\t"<<node[i].len<<endl;
	}
	cout<<endl;
	*/
	int   backtrace[nodeSize];
	short coverage[nodeSize];


	for(int i=0;i<nodeSize;++i) {
		//order[i] = i;
		backtrace[i] = -1;
		coverage[i] =  0;
	}

	//qsort_r(order,nodeSize,sizeof(uint32_t),compare_riAti,&node[0]);


	//node[i].read_seq  < WAITINGLEN + node[j].read_seq + node[j].len;
	//node[j].ref_seq + node[j].len < node[i].ref_seq;

	//connect to head
	for (int i=0;i<nodeSize;++i) {
		if (node[i].read_seq <= waitingLen) {
			backtrace[i] = 0;
			++counterP;
		}
		else
			break;
	}
/*
	if (0 == counterP) {
		for (int i=0;i<nodeSize;++i) {
			if (node[i].read_seq <= WAITINGLENLIMIT) {
				backtrace[i] = 0;
				++counterP;
			}
			else
				break;
		}
	}
*/
	if ( 0 == counterP ) return -1;
	//connect middle part
	for (int i=1; i< nodeSize - 1; i++) {
		//counterP = 0;
		for (int j= i-1;j>=1;--j) {
			if (node[i].read_seq >= node[j].read_seq + node[j].len) {
				if (node[i].read_seq <= waitingLen + node[j].read_seq + node[j].len) {
					if (node[j].ref_seq + node[j].len <= node[i].ref_seq && node[i].ref_seq <= node[j].ref_seq + node[j].len + waitingLen ) {
					//calculate coverage;
						//++counterP;
						if (coverage[i] < coverage[j] + node[j].len) {
							coverage[i] = coverage[j] + node[j].len;
							backtrace[i] = j;
							//cout<<backtrace[i]<<"->"<<i<<"\t"<<coverage[i]<<endl;
						}
					}

				} else 	break;
			}
		}
/*
		if (0==counterP) {
			for (int j= i-1;j>=0;--j) {
				if (node[i].read_seq >= node[j].read_seq + node[j].len) {
					if (node[i].read_seq <= WAITINGLENLIMIT + node[j].read_seq + node[j].len) {
						if (node[j].ref_seq + node[j].len <= node[i].ref_seq) {
							//calculate coverage;
							++counterP;
							if (coverage[i] < coverage[j] + node[j].len) {
								coverage[i] = coverage[j] + node[j].len;
								backtrace[i] = j;
								//cout<<backtrace[i]<<"->"<<i<<"\t"<<coverage[i]<<endl;
							}
						}

					} else 	break;
				}
			}
		}
*/
		//if (0 == counterP) return -1;

	}


	//connnect to tail
	counterP = 0;
	for (int i=nodeSize -2;i>=1;--i) {
		if (node[i].read_seq + waitingLen + node[i].len > node[nodeSize -1].read_seq ) {
			++counterP;
			if (coverage[nodeSize-1] < coverage[i] + node[i].len) {
				coverage[nodeSize-1] = coverage[i] + node[i].len;
				backtrace[nodeSize-1] = i;

			}
		} else break;
	}
/*
	if (0 == counterP) {
		for (int i=nodeSize -2;i>=0;--i) {
			if (node[i].read_seq + WAITINGLENLIMIT + node[i].len > node[nodeSize -1].read_seq) {
				++counterP;
				if (coverage[nodeSize-1] < coverage[i] + node[i].len) {
					coverage[nodeSize-1] = coverage[i] + node[i].len;
					backtrace[nodeSize-1] = i;

				}
			} else break;
		}
	}
*/

	if (0 == counterP) return -1;

	int 			rightOrder[nodeSize];
	int   			counterofRightOrder;

	//cout<<"backtrace: "<<nodeSize-1<<"->";
	rightOrder[0] = nodeSize - 1;
	counterofRightOrder = 1;
	for(int k= backtrace[nodeSize-1];k!=0;k=backtrace[k]) {
		if (k<0) {
			return -1;
			//break;
		} else
			rightOrder[counterofRightOrder++] = k;
			//cout<<k<<"->";
	}

	//for (int k=counterofRightOrder -1;k>=0;--k) {
	//	cout<<rightOrder[k]<<"->";
	//}
	//cout<<endl;
	//for (int i=counterofRightOrder-1;i>=1;--i)

	//fprintf(stderr,"%d",counterofRightOrder);
	//return 1;
	return CalEditDistancewithCigar(rightOrder, counterofRightOrder, read, lenRead, ref, lenRef,left_start,rc,startPos,countChr,sam, countSam, matrix, gapo, gape);

	//get the right order of sequence (rightOrder[] and a counter )
	//ask for revread[200] revref[200] revcigarbuf[401]
	//deal revcigarbuf give its value to cigarbuf, and remember its length,
	//send cigarbuf + len, and remember its usage
	//deal with the first one start recycle
	//first get reflen, get readlen
}
int 	Graphic::dealGraph(uint32_t lenRef, uint32_t lenRead, char *read, char *ref, int *score, uint32_t waitingLen, uint32_t left_start,
		bool rc, uint32_t *startPos, int countChr, SvSam_Rec *sam, int countSam, int8_t *matrix, int gapo, int gape)
{
	vertex start,end;
	start.read_seq = 0;
	start.ref_seq = 0;
	start.len = 0;

	end.read_seq = lenRead - 1;
	end.ref_seq = lenRef - 1;
	end.len = 0;

	int counterP = 0;


	sort(node.begin(),node.end());

	node.insert(node.begin(),start);//maybe changed
	node.push_back(end);

	int nodeSize = node.size();
	//cout<<nodeSize<<endl;
	//prseq(read,lenRead,true);
	//prseq(ref,lenRef,true);

	/*cout<<"After sorting:"<<endl;
	for (int i=0;i<nodeSize;++i) {
		cout<<node[i].read_seq<<"\t"<<node[i].ref_seq<<"\t"<<node[i].len<<endl;
	}
	cout<<endl;
	*/
	int   backtrace[nodeSize];
	short coverage[nodeSize];


	for(int i=0;i<nodeSize;++i) {
		//order[i] = i;
		backtrace[i] = -1;
		coverage[i] =  0;
	}

	//qsort_r(order,nodeSize,sizeof(uint32_t),compare_riAti,&node[0]);


	//node[i].read_seq  < WAITINGLEN + node[j].read_seq + node[j].len;
	//node[j].ref_seq + node[j].len < node[i].ref_seq;

	//connect to head
	for (int i=0;i<nodeSize;++i) {
		if (node[i].read_seq <= waitingLen) {
			backtrace[i] = 0;
			++counterP;
		}
		else
			break;
	}
/*
	if (0 == counterP) {
		for (int i=0;i<nodeSize;++i) {
			if (node[i].read_seq <= WAITINGLENLIMIT) {
				backtrace[i] = 0;
				++counterP;
			}
			else
				break;
		}
	}
*/
	if ( 0 == counterP ) return -1;
	//connect middle part
	for (int i=1; i< nodeSize - 1; i++) {
		//counterP = 0;
		for (int j= i-1;j>=1;--j) {
			if (node[i].read_seq >= node[j].read_seq + node[j].len) {
				if (node[i].read_seq <= waitingLen + node[j].read_seq + node[j].len) {
					if (node[j].ref_seq + node[j].len <= node[i].ref_seq && node[i].ref_seq <= node[j].ref_seq + node[j].len + waitingLen ) {
					//calculate coverage;
						//++counterP;
						if (coverage[i] < coverage[j] + node[j].len) {
							coverage[i] = coverage[j] + node[j].len;
							backtrace[i] = j;
							//cout<<backtrace[i]<<"->"<<i<<"\t"<<coverage[i]<<endl;
						}
					}

				} else 	break;
			}
		}
/*
		if (0==counterP) {
			for (int j= i-1;j>=0;--j) {
				if (node[i].read_seq >= node[j].read_seq + node[j].len) {
					if (node[i].read_seq <= WAITINGLENLIMIT + node[j].read_seq + node[j].len) {
						if (node[j].ref_seq + node[j].len <= node[i].ref_seq) {
							//calculate coverage;
							++counterP;
							if (coverage[i] < coverage[j] + node[j].len) {
								coverage[i] = coverage[j] + node[j].len;
								backtrace[i] = j;
								//cout<<backtrace[i]<<"->"<<i<<"\t"<<coverage[i]<<endl;
							}
						}

					} else 	break;
				}
			}
		}
*/
		//if (0 == counterP) return -1;

	}


	//connnect to tail
	counterP = 0;
	for (int i=nodeSize -2;i>=1;--i) {
		if (node[i].read_seq + waitingLen + node[i].len > node[nodeSize -1].read_seq ) {
			++counterP;
			if (coverage[nodeSize-1] < coverage[i] + node[i].len) {
				coverage[nodeSize-1] = coverage[i] + node[i].len;
				backtrace[nodeSize-1] = i;

			}
		} else break;
	}
/*
	if (0 == counterP) {
		for (int i=nodeSize -2;i>=0;--i) {
			if (node[i].read_seq + WAITINGLENLIMIT + node[i].len > node[nodeSize -1].read_seq) {
				++counterP;
				if (coverage[nodeSize-1] < coverage[i] + node[i].len) {
					coverage[nodeSize-1] = coverage[i] + node[i].len;
					backtrace[nodeSize-1] = i;

				}
			} else break;
		}
	}
*/

	if (0 == counterP) return -1;

	int 			rightOrder[nodeSize];
	int   			counterofRightOrder;

	//cout<<"backtrace: "<<nodeSize-1<<"->";
	rightOrder[0] = nodeSize - 1;
	counterofRightOrder = 1;
	for(int k= backtrace[nodeSize-1];k!=0;k=backtrace[k]) {
		if (k<0) {
			return -1;
			//break;
		} else
			rightOrder[counterofRightOrder++] = k;
			//cout<<k<<"->";
	}

	//for (int k=counterofRightOrder -1;k>=0;--k) {
	//	cout<<rightOrder[k]<<"->";
	//}
	//cout<<endl;
	//for (int i=counterofRightOrder-1;i>=1;--i)

	//fprintf(stderr,"%d",counterofRightOrder);
	//return 1;
	return CalEditDistancewithCigar(rightOrder, counterofRightOrder, read, lenRead, ref, lenRef,left_start,rc,startPos,countChr,sam, countSam, matrix, gapo, gape);

	//get the right order of sequence (rightOrder[] and a counter )
	//ask for revread[200] revref[200] revcigarbuf[401]
	//deal revcigarbuf give its value to cigarbuf, and remember its length,
	//send cigarbuf + len, and remember its usage
	//deal with the first one start recycle
	//first get reflen, get readlen
}
