## DESCRIPTION
# This is generic brick that can be use for all kind of purpose.
# The only requirements are to declare some intput file(s), one output file
# and a command line to be exected to create the output from the input

## VARIABLES
%%_COMMAND_LINE	= # this is the command line that will be run
%%_INPUT_FILE		=# Specify input files if the command needs to be re-run when intput files changes
%%_OUTPUT_FILE  =# Specify an output file if the command produce some file, this is

.PHONY: %%

%%: $(%%_OUTPUT_FILE) $(%%_INPUT_FILE)
ifeq ($(strip $(%%_OUTPUT_FILE)),)
	$(%%_COMMAND_LINE)
endif

%%_clean:
ifneq ($(strip $(%%_OUTPUT_FILE)),)
	-rm $(%%_OUTPUT_FILE)
endif
	
$(%%_OUTPUT_FILE): $(%%_INPUT_FILE)
	$(%%_COMMAND_LINE)
