# INPUT VARS
%%_RAWS = undef
%%_INDEX = undef
%%_INPUT_DIRECTORY = raw
%%_RAW_EXTENSION = fastq
%%_BINARY = crac
%%_VERSION = $(shell $(%%_BINARY) -v | head -1 | cut -d " " -f 3)
%%_OUTPUT_DIRECTORY =	$(core_MAPPING_DIRECTORY)/%%-$(%%_VERSION)
%%_KMER_LENGTH = 22
%%_REMOVE_SAMS = false
%%_OPTIONS =

# OUTPUT VARS
%%_RAWS_FILENAMES = $(notdir $(%%_RAWS))
%%_SAMS = $(%%_RAWS_FILENAMES:%_1.$(%%_RAW_EXTENSION)=$(%%_OUTPUT_DIRECTORY)/%.sam)
%%_BAMS = $(%%_SAMS:.sam=.bam)

.PHONY: %%
%%: $(%%_SAMS)

%%_bam: %% $(%%_SAMS:.sam=.bam)

%%_clean:
	-rm $(%%_SAMS)
	-rm -rf $(%%_OUTPUT_DIRECTORY)

$(%%_OUTPUT_DIRECTORY): $(MAPPING_DIRECTORY)
	mkdir -p $@

$(%%_OUTPUT_DIRECTORY)/%.sam: $(%%_INPUT_DIRECTORY)/%_1.$(%%_RAW_EXTENSION) $(%%_OUTPUT_DIRECTORY)
	$(%%_BINARY) -i $(%%_INDEX) -r $< $(<:_1.$(%%_RAW_EXTENSION)=_2.$(%%_RAW_EXTENSION)) -k $(%%_KMER_LENGTH) $(%%_OPTIONS) -o - --summary $(@:.sam=-summary.txt) > $@  2> $(@:.sam=.log)
