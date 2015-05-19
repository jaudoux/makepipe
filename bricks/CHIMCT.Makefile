## INPUT VARS
%%_SAM_FILE					= undef # A SAM file (or BAM) generated by CRAC
%%_ANNOTATION_FILE	= undef # This is the GFF file that we will use for annotating chimera

## OUTPUT VARS
%%_OUTPUT_FILE			= undef # CHIMCT output in tabulated format (see chimCT -h for more informations)

## OPTIONS
%%_BINARY						= chimCT # Binary used to call CHIMCT
%%_OPTIONS 					= # Options provided to CHIMCT

.PHONY: %%
%%: $(%%_OUTPUT_FILE) 

%%_clean:
	-rm $(basename $(%%_OUTPUT_FILE))*
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

$(%%_OUTPUT_FILE): $(%%_SAM_FILE)
	@mkdir -p $(dir $@)
	$(%%_BINARY) -n $(basename $(notdir $<)) -g $(%%_ANNOTATION_FILE) -s $< $(%%_OPTIONS) > $@ 2> $(@:$(suffix $@)=.log)
