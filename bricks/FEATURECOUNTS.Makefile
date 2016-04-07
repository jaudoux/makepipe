## INPUT VARS
# These are the variables that FEATURE COUNTS needs in input
%%_ANNOTATION_FILE 		= undef# This is the annotation file that we will use to count (in GTF format)
%%_INPUT_ALIGNMENTS		= undef# This is the Alignement file (SAM or BAM)

## OUTPUT VARS
%%_OUTPUT_FILE		 		= undef# This is the count file that will be generated

## OPTIONS
%%_BINARY							= featureCounts# Binary used to run FEATURECOUNTS
%%_OPTIONS						= # Options that can be specified to FEATURECOUNTS (run featurecounts -h for more informations

%%_CLEAN_COUNT_TABLE = true# If this is activated the feature count output is post-processed to only print a nice count-table without other informations

%%_SAMPLE_NAMES = $(notdir $(basename $(%%_INPUT_ALIGNMENTS)))#Name given to the samples if the CLEAN_COUT_TABLE is set to true
%%_OLD_COUNT_TABLE = $(addsuffix .original,$(%%_OUTPUT_FILE))#Filename given to the original count table given by featureCouts if CLEAN_COUNT_TABLE is activated

ifeq ($(%%_CLEAN_COUNT_TABLE),true)
%%_FEATURECOUNTS_OUTPUT_FILE = $(%%_OUTPUT_FILE).original
else
%%_FEATURECOUNTS_OUTPUT_FILE = $(%%_OUTPUT_FILE)
endif

empty:=
space:= $(empty) $(empty)
tab:=\t

%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE) $(%%_OUTPUT_FILE).summary
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

$(%%_FEATURECOUNTS_OUTPUT_FILE): $(%%_INPUT_ALIGNMENTS) $(%%_ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) -a $(%%_ANNOTATION_FILE) -o $@ $(%%_INPUT_ALIGNMENTS)

$(%%_OUTPUT_FILE): $(%%_FEATURECOUNTS_OUTPUT_FILE)
	echo -e '$(subst $(space),$(tab),feature $(%%_SAMPLE_NAMES))' > $@
	tail -n +3 $< | cut -f1,7- >> $@
