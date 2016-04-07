## DESCRIPTION
# This brick execute the "kallisto" software
# (https://pachterlab.github.io/kallisto/) wich is a program for quantifying
# abundances of transcripts from RNA-Seq data.

## INPUT VARS
%%_READS 			= undef#Reads in FASTQ
%%_INDEX 			= undef#Kallisto index
%%_OUTPUT_DIR = undef#Output directory

## OPTIONS
%%_BINARY			= kallisto#The binary used to execute Kallisto
%%_NB_THREADS = 1#Number of threads to be used
%%_OPTIONS		= #Other options to be passed to kallisto

## COMPUTED VARS
%%_ABUNDANCE_FILE			= $(%%_OUTPUT_DIR)/abundance.tsv#
%%_ABUNDANCE_H5_FILE	= $(%%_OUTPUT_DIR)/abundance.h5#
%%_RUN_INFO_FILE			= $(%%_OUTPUT_DIR)/run_info.json#


%%_PRODUCED_FILES			= $(%%_ABUNDANCE_FILE) $(%%_ABUNDANCE_H5_FILE) $(%%_RUN_INFO_FILE)
%%_OPTIONS += -i $(%%_INDEX) -o $(%%_OUTPUT_DIR) -t $(%%_NB_THREADS)

.PHONY: %%

%%: $(%%_ABUNDANCE_FILE)

$(%%_ABUNDANCE_H5_FILE) $(%%_RUN_INFO_FILE): $(%%_ABUNDANCE_FILE)

%%_clean:
	-rm $(%%_PRODUCED_FILES)
	-rmdir -p --ignore-fail-on-non-empty $(%%_OUTPUT_DIR)

#.SERIAL: $(%%_ABUNDANCE_FILE)
$(%%_ABUNDANCE_FILE): $(%%_READS) $(%%_INDEX)
	mkdir -p $(dir $@)
	$(%%_BINARY) quant $(%%_OPTIONS) $(%%_READS) $(%%_PIPE_COMMANDS)
