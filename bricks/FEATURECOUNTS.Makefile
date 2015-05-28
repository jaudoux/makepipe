## INPUT VARS
# These are the variables that FEATURE COUNTS needs in input
%%_ANNOTATION_FILE 		= undef # This is the annotation file that we will use to count (in GTF format)
%%_INPUT_ALIGNMENTS		= undef # This is the Alignement file (SAM or BAM)

## OUTPUT VARS
%%_OUTPUT_FILE		 		= undef # This is the count file that will be generated

## OPTIONS
%%_BINARY							= featureCounts # Binary used to run FEATURECOUNTS
%%_OPTIONS						= # Options that can be specified to FEATURECOUNTS (run featurecounts -h for more informations

%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE) $(%%_OUTPUT_FILE).summary
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

$(%%_OUTPUT_FILE): $(%%_INPUT_ALIGNMENTS) $(%%_ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) -a $(%%_ANNOTATION_FILE) -o $@ $(%%_INPUT_ALIGNMENTS) > $@
