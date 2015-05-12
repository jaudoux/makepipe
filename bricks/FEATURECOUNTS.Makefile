%%_ANNOTATION_FILE 		= undef
%%_INPUT_ALIGNMENTS		= undef
%%_OUTPUT_DIRECTORY 	=	$(core_QUANTIFICATION_DIRECTORY)/%%
%%_OUTPUT_FILENAME 		= gene_counts.tsv
%%_BINARY							= featureCounts
%%_OPTIONS						=

%%_OUTPUT_FILE = $(addprefix $(%%_OUTPUT_DIRECTORY)/, $(%%_OUTPUT_FILENAME))

%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE)
	-rm -rf $(%%_OUTPUT_DIRECTORY)

$(%%_OUTPUT_FILE): $(%%_INPUT_ALIGNMENTS) $(%%_ANNOTATION_FILE)
	@mkdir -p $(%%_OUTPUT_DIRECTORY)
	$(%%_BINARY) $(%%_OPTIONS) -a $(%%_ANNOTATION_FILE) -o $@ $(%%_INPUT_ALIGNMENTS) > $@
