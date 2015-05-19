# INPUT VARS
%%_INDEX 							= undef
%%_READS 							= undef

# OUTPUT VARS
%%_OUTPUT_FILE				= undef

# OPTIONS
%%_BINARY 						= crac
%%_KMER_LENGTH 				= 22
%%_OPTIONS 						=

# COMPUTED VARS
%%_VERSION 						= $(shell $(%%_BINARY) -v | head -1 | cut -d " " -f 3)
%%_OPTIONS += -i $(%%_INDEX) -k $(%%_KMER_LENGTH) --bam
%%_OUTPUT_EXTENSION 	= $(suffix $(%%_OUTPUT_FILE))

#%%_RAWS = undef
#%%_READS_DIRECTORY 		= raw
#%%_READS_EXTENSION 		= .fastq
#%%_FIRST_PAIR_SUFFIX 	= _1
#%%_SECOND_PAIR_SUFFIX	= _2
#%%_OUTPUT_DIRECTORY 	=	$(core_MAPPING_DIRECTORY)/%%-$(%%_VERSION)
#%%_OUTPUT_EXTENSION 	= .bam
#%%_REMOVE_SAMS 				= false

# OUTPUT VARS
#%%_RAWS_FILENAMES = $(notdir $(%%_RAWS))
#%%_SAMS 							= $(addprefix $(%%_OUTPUT_DIRECTORY)/, $(addsuffix $(%%_OUTPUT_EXTENSION), $(%%_SAMPLES)))

# COMPUTED VARS

.PHONY: %%
%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE)
	#-rm -rf $(%%_OUTPUT_DIRECTORY)

.SERIAL: $(%%_OUTPUT_FILE)
$(%%_OUTPUT_FILE): $(%%_READS)
	mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) -r $(wordlist 1, 2, $^) -o - --summary $(@:$(%%_OUTPUT_EXTENSION)=-summary.txt) 2> $(@:$(%%_OUTPUT_EXTENSION)=.log) | samtools sort - $(@:%$(%%_OUTPUT_EXTENSION)=%)
	samtools index $@
