## DESCRIPTION
# HISAT is a Mapping softwares that is able to align spliced reads to a reference genome.

## SYNOPSIS
# In order to run HISAT, one only need to specify a HISAT index, sequence files 
# (in FASTA/FASTQ) and a prefix for the output files.

## INPUT VARS:
%%_INDEX 							= undef# Filename of the index to be used
%%_READS 							= undef# List of reads files (FASTA or FASTQ) to be passed to HISAT

## OUTPUT VARS
%%_OUTPUT_PREFIX			= hisat# Prefix for name of the output files (see below)
%%_OUTPUT_FILE				= $(%%_OUTPUT_PREFIX).bam# SAM/BAM file generated by HISAT
%%_LOG_FILE       		= $(%%_OUTPUT_PREFIX).log# Log file from STDERR
%%_NOVEL_SPLICE_FILE	= # Output file for novel splice sites

## OPTIONS
%%_BINARY 						= hisat# Binary to use for running HISAT
%%_NB_THREADS					= 1# Number of threads allocated to the cractools binary
%%_OPTIONS 						= # Options that can be specified to HISAT (run `hisat` for more information) 
%%_SORT_BAM						= true# If true, output bam will be piped to samtools for sorting
%%_INDEX_BAM					= true# If true, sorted bam will be indexed by samtools to create a .bam.bai index file
%%_SAMTOOLS_OPTIONS   = -@4# Options passed to samtools view

## COMPUTED VARS
%%_VERSION 						= $(shell $(%%_BINARY) --version 2>&1 | head -1 | cut -f3 -d " ") # HISAT version

%%_OPTIONS += -x $(%%_INDEX) -p $(%%_NB_THREADS)# Crac Options
%%_PRODUCED_FILES = $(%%_OUTPUT_FILE) $(%%_NOVEL_SPLICE_FILE) $(%%_LOG_FILE)

ifeq ($(words $(%%_READS)),2)
%%_OPTIONS += -1 $(word 1,$(%%_READS)) -2 $(word 2,$(%%_READS))
else
%%_OPTIONS += -U $(%%_READS)
endif

ifneq ($(%%_NOVEL_SPLICE_FILE),)
%%_OPTIONS += --novel-splicesite-outfile $(%%_NOVEL_SPLICE_FILE)
endif

%%_PIPE_COMMANDS = > $(%%_OUTPUT_FILE)
ifeq ($(%%_SORT_BAM),true)
%%_PIPE_COMMANDS = | samtools view -bS - | samtools sort $(%%_SAMTOOLS_OPTIONS) - $(basename $(%%_OUTPUT_FILE))
ifeq ($(%%_INDEX_BAM),true)
%%_PIPE_COMMANDS += && samtools index $(%%_OUTPUT_FILE)
%%_PRODUCED_FILES += $(basename $(%%_OUTPUT_FILE)).bam.bai
endif
endif

.PHONY: %%

%%: $(%%_OUTPUT_FILE)

$(%%_NOVEL_SPLICE_FILE): $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_PRODUCED_FILES)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

.SERIAL: $(%%_OUTPUT_FILE)
$(%%_OUTPUT_FILE): $(%%_READS) $(%%_INDEX).1.ht2
	mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) 2> $(%%_LOG_FILE) $(%%_PIPE_COMMANDS)
