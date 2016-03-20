## DESCRIPTION
# This brick uses Jellyfish to count k-mers in one or more FAST[A,Q] files and
# save the k-mer database on disk.

## INPUT VARS
%%_READS				=	undef#Reads file (FASTQ) to index
%%_KMER_LENGTH	=	30#K-mer length to use

## OUTPUT VARS
%%_OUTPUT_FILE	= mer_counts.jf#Output produced bu Jellyfish

## OPTIONS
%%_BINARY 						= jellyfish# Binary to use for running Jellyfish
%%_INITIAL_HASH_SIZE 	= 1000#Initial size of jellyfish hash 
%%_LOWER_COUNT				= 2#Do not output k-mers with counts lower than this threshold
%%_NB_THREADS					= 1#Number of threads to use
%%_OPTIONS						= -C#Options given to jellyfish

%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

.SERIAL: $(%%_OUTPUT_FILE)
$(%%_OUTPUT_FILE): $(%%_READS)
	mkdir -p $(dir $@)
	$(%%_BINARY) count -t $(%%_NB_THREADS) -m $(%%_KMER_LENGTH) -L $(%%_LOWER_COUNT) -s $(%%_INITIAL_HASH_SIZE) -o $@ $(%%_OPTIONS) <(zcat $^)

