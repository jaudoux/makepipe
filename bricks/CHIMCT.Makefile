# INPUT VARS
%%_SAM_FILE					= undef
%%_ANNOTATION_FILE	= undef

# OUTPUT VARS
%%_OUTPUT_FILE			= undef

# OPTIONS
%%_BINARY						= chimCT
%%_OPTIONS 					=

#%%_SAM_EXTENSION		= bam
#%%_OUTPUT_DIRECTORY	=	$(core_CHIMERAS_DIRECTORY)/%%
#%%_SAMS_FILENAMES 	= $(basename $(notdir $(%%_SAMS)))
#%%_CSVS_EXTENSION		= -chimeras.csv
#%%_CSVS 						= $(%%_SAMS_FILENAMES:%=$(%%_OUTPUT_DIRECTORY)/%$(%%_CSVS_EXTENSION))

.PHONY: %%
%%: $(%%_OUTPUT_FILE) 

%%_clean:
	-rm $(%%_OUTPUT_FILE)

$(%%_OUTPUT_FILE): $(%%_SAM_FILE)
	@mkdir -p $(dir $@)
	$(%%_BINARY) -n $(basename $(notdir $<)) -g $(%%_ANNOTATION_FILE) -s $< $(%%_OPTIONS) > $@ 2> $(@:$(suffix $@)=.log)
