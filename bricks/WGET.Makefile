## DESCRIPTION
# This brick download a file from a distant location given an URL (http/ftp).
# The WGET program is used for that purpose (see C<man wget>).

## INPUT VARS
%%_URL					= undef#URL used for the download

## OUTPUT VARS
%%_OUTPUT_FILE  = undef#Name of the output file

%%: $(%%_OUTPUT_FILE)

%%_clean:
	rm $(%%_OUTPUT_FILE)

$(%%_OUTPUT_FILE):
	mkdir -p $(dir $@)
	wget -O $@ $(%%_URL)
