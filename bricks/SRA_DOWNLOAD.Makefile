## DESCRIPTION
# This brick downloads files from SRA using an SRR ids by using the 'fastq-dump' software
# from SRA toolkit.

## INPUT VARS
%%_SRR_ID							=	undef# SRA identifiers for the run (with SRR prefix)
%%_TYPE   						= PE#Type of sequencing : PE (Paired-End) or SE (Single-End)

## OPTIONS
%%_MAX_READS 					=	#Maximum number of reads to download
%%_GZIP_OUTPUT				= true#Output is gzipped (true or false)
%%_FASTQ_DUMP_OPTIONS	=	#Other options to give to fastq-dump

## OUTPUT VARS
%%_READS_PREFIX				=	$(%%_SRR_ID)

## COMPUTED VARS
%%_READS							=	#Reads files downloaded

# Manage READS output
ifeq ($(%%_TYPE),PE) 
%%_READS							=	$(%%_READS_PREFIX)_1.fastq.gz $(%%_READS_PREFIX)_2.fastq.gz
%%_FASTQ_DUMP_OPTIONS += --split-files
else
%%_READS							=	$(%%_READS_PREFIX).fastq.gz
endif

ifeq ($(%%_GZIP_OUTPUT),true)
%%_FASTQ_DUMP_OPTIONS += --gzip
endif

ifneq ($(%%_MAX_READS),)
%%_FASTQ_DUMP_OPTIONS += -X $(%%_MAX_READS)
endif

.PHONY: %%
%%: $(%%_READS)

%%_clean:
	rm $(%%_READS)

$(word 1, $(%%_READS)):
	fastq-dump $(%%_FASTQ_DUMP_OPTIONS) $(%%_SRR_ID)
	mkdir -p $(dir $(%%_READS_PREFIX))
	rename 's/$(%%_SRR_ID)/$(subst /,\/,$(%%_READS_PREFIX))/' $(%%_SRR_ID)*

$(word 2, $(%%_READS)): $(word 1, $(%%_READS))
