## INPUT VARS
# These are the variables that FEATURE COUNTS needs in input
%%_ANNOTATION_FILE 		= undef # This is the annotation file that we will use to count (in GTF format)
%%_INPUT_ALIGNMENTS		= undef # This is the Alignement file (SAM or BAM)

## OUTPUT VARS
%%_OUTPUT_FILE		 		= undef # This is the count file that will be generated

## OPTIONS
%%_BINARY							= featureCounts # Binary used to run FEATURECOUNTS
%%_OPTIONS						= # Options that can be specified to FEATURECOUNTS (run featurecounts -h for more informations

#COMPUTED VARS

#%%_OUTPUT_FILENAME 		= gene_counts.tsv
#%%_OUTPUT_DIRECTORY 	=	$(core_QUANTIFICATION_DIRECTORY)/%%
#%%_OUTPUT_FILE = $(addprefix $(%%_OUTPUT_DIRECTORY)/, $(%%_OUTPUT_FILENAME))

%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE)
	-rm -rf $(%%_OUTPUT_DIRECTORY)

$(%%_OUTPUT_FILE): $(%%_INPUT_ALIGNMENTS) $(%%_ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) -a $(%%_ANNOTATION_FILE) -o $@ $(%%_INPUT_ALIGNMENTS) > $@
