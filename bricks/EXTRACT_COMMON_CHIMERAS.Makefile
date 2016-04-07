## INPUT VARS
%%_CHIMCT_FILE			= undef # This is the CHIMCT output file(s) 

## OUTPUT VARS
%%_OUTPUT_FILE			= undef # This is the results file that will be generated

## OPTIONS
%%_BINARY	= extractCommonChimeras.pl # Binary used to call EXTRACT_COMMON_CHIMERAS
%%_OPTIONS	= # Options provided to EXTRACT_COMMON_CHIMERAS
%%_SUMMARY_OUTPUT 	=  # Output for making some statistics
%%_RDATA_OUTPUT		=  # Output to format the statistics for R graphics
%%_SUMMARY = # If correspondant OUTPUT not empty, make some statitics
%%_RDATA = # If correspondant OUTPUT not empty, format the statistics for R graphics

ifneq ($(%%_SUMMARY_OUTPUT),)
    %%_SUMMARY += --summary
endif
ifneq ($(%%_RDATA_OUTPUT),)
    %%_RDATA += --Rdata
endif

.PHONY: %%
%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE) $(%%_SUMMARY_OUTPUT) $(%%_RDATA_OUTPUT)

$(%%_OUTPUT_FILE): $(%%_CHIMCT_FILE)
	$(%%_BINARY) $(%%_OPTIONS) $(%%_SUMMARY) $(%%_SUMMARY_OUTPUT) $(%%_RDATA) $(%%_RDATA_OUTPUT) $^ > $@
