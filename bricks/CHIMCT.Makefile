# CHIMCT
%%_SAMS							= undef
%%_INPUT_DIRECTORY  = undef
%%_GFF							= undef
%%_SAM_EXTENSION		= bam
%%_BINARY						= chimCT
%%_OUTPUT_DIRECTORY	=	$(core_CHIMERAS_DIRECTORY)/%%
%%_SAMS_FILENAMES 	= $(basename $(notdir $(%%_SAMS)))
%%_CSVS_EXTENSION		= -chimeras.csv
%%_CSVS 						= $(%%_SAMS_FILENAMES:%=$(%%_OUTPUT_DIRECTORY)/%$(%%_CSVS_EXTENSION))
%%_OPTIONS 					=

.PHONY: %%
%%: $(%%_CSVS) 

%%_clean:
	-rm $(%%_CSVS)
	-rm -rf $(%%_OUTPUT_DIRECTORY) 

$(%%_OUTPUT_DIRECTORY)/%$(%%_CSVS_EXTENSION): $(%%_INPUT_DIRECTORY)/%.$(%%_SAM_EXTENSION) 
	@mkdir -p $(%%_OUTPUT_DIRECTORY)
	$(%%_BINARY) -n $(basename $(notdir $<)) -g $(%%_GFF) -s $< $(%%_OPTIONS) > $@ 2> $(@:.csv=.log)
