#ifndef HASH_H_
#define HASH_H_
	#define NUM_FILE 256
	#define N_LIMIT 9
	#include <stdint.h>
	typedef struct{
		uint32_t *pointer;
		//char 	**seq_bkt;
		uint32_t *seq_bkt;
	}Hashtab;

	typedef struct{
		uint32_t sed_value;
		uint32_t sed_bkt_pos;
	}seed;

	class Hash{	
	public:
		//Hashtab *gen_hash(char *genome,uint32_t len_genome, char *path, uint32_t len_sed);
		Hashtab *load_hashfile(char *path, uint32_t len_genome, uint32_t len_sed);
		Hashtab *write_hashfile(char *path,char *genome,uint32_t len_genome,uint32_t len_sed);
	};
	extern uint8_t trans[];
	int transfer(char *genome,uint32_t *l2r,uint32_t len_sed);
#endif
