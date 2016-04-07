## DESCRIPTION

## VARIABLES
%%_COMMAND_LINE	= # this is the command line that will be run
%%_INPUT_PATTERN	=undef# Specify input files if the command needs to be re-run when intput files changes
%%_OUTPUT_PATTERN  =undef# Specify an output file if the command produce some file, this is

%$(%%_OUTPUT_PATTERN): %$(%%_INPUT_PATTERN)
	$(%%_COMMAND_LINE)
