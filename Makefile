CRAC_INDEX = /data/indexes/crac/GRCh37
CRAC_OPTIONS = --stringent-chimera
CRAC_KMER_LENGTH = 22
CRAC_INPUT_DIRECTORY = raw
CRAC_RAWS = raw/test1_1.fastq raw/test2_1.fastq
CRAC_BINARY = /data/bin/crac

CHIMCT_SAMS = $(CRAC_SAMS)
CHIMCT_INPUT_DIRECTORY = $(CRAC_OUTPUT_DIRECTORY)

all: crac_bam chimCT

#chimCT
#$(CRAC_SAMS:.sam=.bam)

clean :
	rm $(CRAC_SAMS) $(CHIMCT_CSVS)

include Makefile.global
