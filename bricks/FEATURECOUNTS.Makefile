# INPUT VARS
%%_ANNOTATION_FILE 		= undef
%%_INPUT_ALIGNMENTS		= undef

# OUTPUT VARS
%%_OUTPUT_FILE		 		= undef

# OPTIONS
%%_BINARY							= featureCounts
%%_OPTIONS						=

# COMPUTED VARS

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
