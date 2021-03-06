## DESCRIPTION
# CRAC is a Mapping softwares that is able to align spliced reads to a reference genome.

## SYNOPSIS
# In order to run CRAC, one only need to specify a CRAC index, sequence files 
# (in FASTA/FASTQ) and a prefix for the output files.

## INPUT VARS:
%%_INDEX 							= undef# Filename of the index to be used
%%_READS 							= undef# List of reads files (FASTA or FASTQ) to be passed to CRAC

## OUTPUT VARS
%%_OUTPUT_PREFIX			= crac# Prefix for name of the output files (see below)
%%_OUTPUT_FILE				= $(%%_OUTPUT_PREFIX).bam# SAM/BAM file generated by CRAC
%%_SUMMARY_FILE       = $(%%_OUTPUT_PREFIX)-summary.txt# Summary file generated by CRAC (--summary option)	
%%_LOG_FILE       		= $(%%_OUTPUT_PREFIX).log# Log file from STDERR

## OPTIONS
%%_BINARY 						= crac# Binary to use for running CRAC
%%_KMER_LENGTH 				= 22# Size of the k-mer used by CRAC
%%_NB_THREADS					= 1# Number of threads allocated to the cractools binary
%%_OPTIONS 						= # Options that can be specified to CRAC (run `crac -f` for more information) 
%%_SORT_BAM						= true# If true, output bam will be piped to samtools for sorting
%%_INDEX_BAM					= true# If true, sorted bam will be indexed by samtools to create a .bam.bai index file
%%_SAMTOOLS_OPTIONS   = -@4# Options passed to samtools view

## COMPUTED VARS
%%_VERSION 						= $(shell $(%%_BINARY) -v | head -1 | cut -d " " -f 3)# Crac version
%%_OPTIONS += -i $(%%_INDEX) -k $(%%_KMER_LENGTH) --bam --nb-tags-info-stored 10000# Crac Options
%%_PRODUCED_FILES = $(%%_OUTPUT_FILE) $(%%_SUMMARY_FILE) $(%%_LOG_FILE)


%%_PIPE_COMMANDS = > $(%%_OUTPUT_FILE)
ifeq ($(%%_SORT_BAM),true)
%%_PIPE_COMMANDS = | samtools sort $(%%_SAMTOOLS_OPTIONS) - -o $(%%_OUTPUT_FILE)
ifeq ($(%%_INDEX_BAM),true)
%%_PIPE_COMMANDS += && samtools index $(%%_OUTPUT_FILE)
%%_PRODUCED_FILES += $(basename $(%%_OUTPUT_FILE)).bam.bai
endif
endif

.PHONY: %%
#%%: $(%%_SUMMARY_FILE)
%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_PRODUCED_FILES)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

$(%%_LOG_FILE) $(%%_SUMMARY_FILE): $(%%_OUTPUT_FILE)

.SERIAL: $(%%_OUTPUT_FILE)
$(%%_OUTPUT_FILE): $(%%_READS)
	mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) --nb-threads $(%%_NB_THREADS) -r $(wordlist 1, 2, $^) -o - --summary $(%%_SUMMARY_FILE) 2> $(%%_LOG_FILE) $(%%_PIPE_COMMANDS)
