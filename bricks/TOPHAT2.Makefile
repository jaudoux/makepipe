## INPUT VARS
%%_INDEX 							= undef# Filename of the index to be used
%%_READS 							= undef# List of reads files (FASTA or FASTQ) to be passed to CRAC

## OUTPUT VARS
%%_OUTPUT_DIR					= tophat2
%%_LOG_FILE						= $(%%_OUTPUT_DIR)/log.txt

## OPTIONS
%%_BINARY							= tophat2
%%_NB_THREADS					= 1
%%_OPTIONS						= 

## AUTOMATED OUTPUT VARS
%%_OUTPUT_FILE				= $(%%_OUTPUT_DIR)/accepted_hits.bam
%%_JUNCTIONS_BED			= $(%%_OUTPUT_DIR)/junctions.bed
%%_INSERTIONS_BED			= $(%%_OUTPUT_DIR)/insertions.bed
%%_DELETIONS_BED			= $(%%_OUTPUT_DIR)/deletions.bed

## COMPUTED VARS
%%_VERSION 						= $(shell $(%%_BINARY) -v | cut -d" " -f2)
%%_OPTIONS += -o $(%%_OUTPUT_DIR) -p $(%%_NB_THREADS)
%%_PRODUCED_FILES = $(%%_OUTPUT_FILE) $(%%_JUNCTIONS_BED) $(%%_INSERTIONS_BED) $(%%_DELETIONS_BED)

.PHONY: %%
%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_PRODUCED_FILES)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

$(%%_JUNCTIONS_BED) $(%%_INSERTIONS_BED) $(%%_DELETIONS_BED): $(%%_OUTPUT_FILE)

.SERIAL: $(%%_OUTPUT_FILE)
$(%%_OUTPUT_FILE): $(%%_READS)
	mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) $(%%_INDEX) $(wordlist 1, 2, $^) 2> $(%%_LOG_FILE) $(%%_PIPE_COMMANDS)
