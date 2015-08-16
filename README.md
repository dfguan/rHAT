###Introduction

rHAT is a seed-and-extension-based noisy long read alignment tool. It is suitable for aligning 3rd generation sequencing reads which are in large read length with relatively high error rate, especially Pacbio's Single Molecule Read-time (SMRT) sequencing reads. 

rHAT indexes the genome with a hash table-based index (regional hash table, RHT) which describes the short tokens occurring in local windows of reference genome. With this index, rHAT adopts a specifically designed seed-and-extension strategy. In the seeding phase, the occurrences of short token matches between partial read and local genomic windows are efficiently calculated to find highly possible sites as candidates for extension. In the extension phase, a sparse dynamic programming-based heuristic approach is adopted for reducing the cost of the alignment between the long noisy read and the local reference sequence. 

rHAT has outstanding throughput on aligning SMRT reads from various prokaryote and eukaryote genomes. Benchmarking on a series of model organism genomes, e.g., E. coli, S. cerevisiae, D. melanogaster, A. thaliana, H. sapiens, etc., demonstrated that it can be two to several times as fast as currently state-of-the-art aligners. Meanwhile, rHAT can sensitively and consecutively aligns the read, i.e., most of the noisy long reads can be end-to-end aligned, and all the bases can be covered.
rHAT is open source and free for non-commercial use.

rHAT is mainly designed by Bo Liu and developed by Dengfeng Guan in Center for Bioinformatics, Harbin Institute of Technology, China.

---

###Memory requirement

The memory usage of rHAT can fit the configurations of most modern servers and workstations. Its peak memory footprint depends on the length of reference genome and the k-mer size parameter setting. We investigated its memory usage for aligning two real SMRT datasets, respectively sequenced from H. sapiens and D. melanogaster genomes, on a server with Intel E5640 CPU at 2.67 GHz, 24 Gigabytes RAM running Linux Ubuntu 10.04. The read were respectively aligned to GRCh37/hg19 and DM5 reference genomes, and the peak memory footprint is as following.

H. sapiens dataset:
```
k-mer size=13 (default setting), 13.70 Gigabytes 
k-mer size=15 (max k-mer size),  17.48 Gigabytes 
```

D. melanogaster dataset:
```
k-mer size=13 (default setting), 1.01  Gigabytes 
k-mer size=15 (max k-mer size),  4.76 Gigabytes 
```

---

###Installation

Current version of rHAT needs to be run on Linux operating system.

The source code is written in C++, and can be directly download from: https://github.com/derekguan/rHAT

The makefile is attached. Use the make command for generating the executable file.

---

###Synopsis
```
rHAT-indexer [-k kmerSize] <HashIndexDir> <Reference>

Index reference in RHT format

rHAT-aligner [-w windowsHits] [-m candidates] [-k kmerSize] [-a match] [-b mismatch]
[-q gapOpen] [-r gapExtension] [-t threadNumber] <HashIndexDir> <ReadFile> <Reference>

Align read to its primitive location in Reference
```

---

###Parameters (could be updated in the future for adding new functions)
```
rHAT-indexer:
-k, --kmer-size      INT    The size of the k-mers extracted from the reference genome for indexing[13].

rHAT-aligner:
-w, --window-hits    INT    The max allowed number of windows hit by a k-mer, if a k-mer
                            hits more than -w genomic windows, rHAT would consider the k-mer 
                            is too repetitive, and discard the k-mer. (default = 1000)

-m, --candidates     INT    The number of candidates for extension, this is one of the
                            major heuristic parameters in rHAT. Setting -m high will let
                            rHAT aligns the read to many local sites, which could affect
                            the throughput, while setting -m low may make too few candidates
                            which could affect the sensitivity and accuracy of the alignment. 
                            Based on the benchmarking on a series of simulated and real datasets 
                            from various genomes, we suggest that setting the -m parameter to 
                            5-10 could reconcile the effectiveness and efficiency in most cases. (default = 5)

-k, --kmer-size      INT    The size of the k-mer extracted from the read for generating short 
                            token matches, note that it needs to be same to the setting of -k parameter
                            in rHAT-indexer. It is not allowed to set -k parameter >15 in current version of rHAT,
                            due to the large usage of RAM. It is also worth noting that, for a large 
                            and repetitive reference genome (e.g., mammalian genomes), setting -k too small,
                            e.g., <11, may affect the alignment, since some k-mers may hit too many genomic 
                            windows and ignored by rHAT according to the limit on the windows hits, 
                            i.e., the -w setting. 

-a, --match          INT    score of match for the alignments in extension phase [1]

-b, --mismatch       INT    mismatch penalty for the alignments in extension phase [5]

-q, --gap-open       INT    gap open penalty for the alignments in extension phase [2]

-r, --gap-extension  INT    gap extension penalty for the alignments in extension phase [1]

-l, --local-kmer     INT    The minimum length of the local matches used for SDP, in the extension phase, only the
                            local mathces longer than -l bp will be utilized for building skeleton of alignment [11]

-t, --threads        INT    number of threads [1]
```

---

###Quick start
```
Genome indexing:

rHAT-indexer Genome_Index_Dir Genome_File

Read alignment:

rHAT-aligner Genome_Index_Dir Fastq_File Genome_File > Sam_File
```

---

###Simulation benchmarking:

We simulated a series of datasets from various genomes, i.e., Escherichia coli (complete genome of the 536 strain), Saccharomyces cerevisiae (build sacCer3), Drosophila melanogaster (build DM3), Arabidopsis thaliana (build TAIR10) and Homo sapiens (build GCRh37/hg19), through PBSim (Ono et. al., 2013, Verison 1.0.3). The read length is analogous to Pacbio P5/C3 release. The average read length is 8000 basepairs. These datasets helped us to evaluate the performance of rHAT.
The datasets have been uploaded to Google Drive, and can be downloaded through the following link:
https://drive.google.com/folderview?id=0Bwibkj8plEJrZFlNOG1rd3hBRWM&usp=sharing

---

###Reference

rHAT: Fast aligning noisy long read with regional hashing. *Under Review*.

---

###Contact

For advising, bug reporting and requiring help, please contact ydwang@hit.edu.cn
