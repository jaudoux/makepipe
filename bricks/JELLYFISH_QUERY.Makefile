## DESCRIPTION
# This brick uses Jellyfish library of counted k-mers (see JELLYFISH_COUNT brick) and a list of k-mer to be queried. For each of the queried k-mer, the count in the library is outputed.

## INPUT VARS
%%_JELLYFISH_DATABASE			=	undef#Reads file (FASTQ) to index
%%_SEQUENCE_FILE					=	unedf#K-mer length to use

## OUTPUT VARS
%%_OUTPUT_FILE	= undef#Output produced by Jellyfish query

## OPTIONS
%%_BINARY 						= jellyfish# Binary to use for running Jellyfish
%%_VERSION						= $(shell $(%%_BINARY) --version | cut -d " " -f 2)
%%OPTIONS							= #Options given to jellyfish

%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

.SERIAL: $(%%_OUTPUT_FILE)
$(%%_OUTPUT_FILE): $(%%_JELLYFISH_DATABASE) $(%%_SEQUENCE_FILE)
	mkdir -p $(dir $@)
	$(%%_BINARY) query $(%%_OPTIONS) -o $(%%_OUTPUT_FILE) -s $(%%_SEQUENCE_FILE) $<

