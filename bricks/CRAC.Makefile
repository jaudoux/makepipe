# INPUT VARS
#%%_RAWS = undef
%%_INDEX 							= undef
%%_SAMPLES 						= undef
%%_READS_DIRECTORY 		= raw
%%_READS_EXTENSION 		= .fastq
%%_PAIRED_END_READS  	= 0
%%_FIRST_PAIR_SUFFIX 	= _1
%%_SECOND_PAIR_SUFFIX	= _2
%%_BINARY 						= crac
%%_VERSION 						= $(shell $(%%_BINARY) -v | head -1 | cut -d " " -f 3)
%%_OUTPUT_DIRECTORY 	=	$(core_MAPPING_DIRECTORY)/%%-$(%%_VERSION)
%%_OUTPUT_EXTENSION 	= .bam
%%_KMER_LENGTH 				= 22
%%_REMOVE_SAMS 				= false
%%_OPTIONS 						=

# OUTPUT VARS
#%%_RAWS_FILENAMES = $(notdir $(%%_RAWS))
%%_SAMS 							= $(addprefix $(%%_OUTPUT_DIRECTORY)/, $(addsuffix $(%%_OUTPUT_EXTENSION), $(%%_SAMPLES)))

# COMPUTED VARS
%%_OPTIONS += -i $(%%_INDEX) -k $(%%_KMER_LENGTH) --bam

.PHONY: %%
%%: $(%%_SAMS)

%%_clean:
	-rm $(%%_SAMS)
	-rm -rf $(%%_OUTPUT_DIRECTORY)

$(%%_OUTPUT_DIRECTORY): $(MAPPING_DIRECTORY)
	mkdir -p $@

ifeq ($(%%_PAIRED_END_READS),1)
$(%%_OUTPUT_DIRECTORY)/%$(%%_OUTPUT_EXTENSION): $(%%_READS_DIRECTORY)/%$(%%_FIRST_PAIR_SUFFIX)$(%%_READS_EXTENSION) $(%%_READS_DIRECTORY)/%$(%%_SECOND_PAIR_SUFFIX)$(%%_READS_EXTENSION)
	@mkdir -p $(%%_OUTPUT_DIRECTORY)
	$(%%_BINARY) $(%%_OPTIONS) -r $(wordlist 1, 2, $^) -o - --summary $(@:$(%%_OUTPUT_EXTENSION)=-summary.txt) 2> $(@:$(%%_OUTPUT_EXTENSION)=.log) | samtools sort - $(@:%$(%%_OUTPUT_EXTENSION)=%)
	samtools index $@
else
$(%%_OUTPUT_DIRECTORY)/%$(%%_OUTPUT_EXTENSION): $(%%_READS_DIRECTORY)/%$(%%_READS_EXTENSION)
	@mkdir -p $(%%_OUTPUT_DIRECTORY)
	$(%%_BINARY) $(%%_OPTIONS) -r $< -o - --summary $(@:$(%%_OUTPUT_EXTENSION)=-summary.txt) 2> $(@:$(%%_OUTPUT_EXTENSION)=.log) | samtools sort - $(@:%$(%%_OUTPUT_EXTENSION)=%)
	samtools index $@
endif
