%%_READS = undef
%%_INDEX = undef
%%_OUTPUT_DIR = undef
#%%_NB_THREADS = 1

## OPTIONS
%%_BINARY=kallisto
%%_OPTIONS=

## COMPUTED VARS
%%_ABUNDANCE_FILE=$(%%_OUTPUT_DIR)/abundance.tsv
%%_ABUNDANCE_H5_FILE=$(%%_OUTPUT_DIR)/abundance.h5
%%_RUN_INFO_FILE=$(%%_OUTPUT_DIR)/run_info.json
%%_PRODUCED_FILES=$(%%_ABUNDANCE_FILE) $(%%_ABUNDANCE_H5_FILE) $(%%_RUN_INFO_FILE)
%%_OPTIONS += -i $(%%_INDEX) -o $(%%_OUTPUT_DIR) -t $(%%_NB_THREADS)

.PHONY: %%

%%: $(%%_ABUNDANCE_FILE)

$(%%_ABUNDANCE_H5_FILE) $(%%_RUN_INFO_FILE): $(%%_ABUNDANCE_FILE)

%%_clean:
	-rm $(%%_PRODUCED_FILES)
	-rmdir -p --ignore-fail-on-non-empty $(%%_OUTPUT_DIR)

#.SERIAL: $(%%_ABUNDANCE_FILE)
$(%%_ABUNDANCE_FILE): $(%%_READS)
	$(%%_BINARY) quant $(%%_OPTIONS) $^ $(%%_PIPE_COMMANDS)
