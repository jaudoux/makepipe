# INPUT VARS
%%_SAM_FILE										= undef

# OUTPUT VARS
%%_CHIMERA_OUTPUT							= undef
%%_SPLICE_OUTPUT							= undef
%%_MUTATION_OUTPUT						= undef

# OPTIONS
%%_OPTIONS										= --coverless-splices
%%_BINARY											= cractools-extract

# COMPUTED VARS
#%%_SPLICES_OUTPUT_DIRECTORY 	= $(core_SPLICES_DIRECTORY)
#%%_MUTATIONS_OUTPUT_DIRECTORY = $(core_MUTATIONS_DIRECTORY)
#%%_CHIMERAS_OUTPUT_DIRECTORY 	= $(core_CHIMERAS_DIRECTORY)
#%%_SPLICES_EXTENSION 					= -splices.bed
#%%_MUTATIONS_EXTENSION				= .vcf
#%%_CHIMERAS_EXTENSION					= -chimeras.tsv
#%%_SAMPLE_NAME								= $(basename $(notdir $(%%_SAM_FILE)))
#%%_CHIMERA_OUTPUT							= $(addprefix $(%%_CHIMERAS_OUTPUT_DIRECTORY)/, $(addsuffix $(%%_CHIMERAS_EXTENSION), $(%%_SAMPLE_NAME)))
#%%_SPLICE_OUTPUT							= $(addprefix $(%%_SPLICES_OUTPUT_DIRECTORY)/, $(addsuffix $(%%_SPLICES_EXTENSION), $(%%_SAMPLE_NAME)))
#%%_MUTATION_OUTPUT							= $(addprefix $(%%_MUTATIONS_OUTPUT_DIRECTORY)/, $(addsuffix $(%%_MUTATIONS_EXTENSION), $(%%_SAMPLE_NAME)))

%%: $(%%_SPLICE_OUTPUT)

%%_clean:
	-rm $(%%_MUTATION_OUTPUT)
	-rm $(%%_SPLICE_OUTPUT)
	-rm $(%%_CHIMERA_OUTPUT)

$(%%_MUTATION_OUTPUT): $(%%_SPLICE_OUTPUT)

$(%%_CHIMERA_OUTPUT): $(%%_SPLICE_OUTPUT)

$(%%_SPLICE_OUTPUT): $(%%_SAM_FILE)
	@mkdir -p $(dir $(%%_CHIMERA_OUTPUT) $(%%_SPLICE_OUTPUT) $(%%_MUTATION_OUTPUT))
	$(%%_BINARY) $(%%_OPTIONS) $< -s $@ -m $(%%_MUTATION_OUTPUT) -c $(%%_CHIMERA_OUTPUT)
	@touch $(%%_MUTATION_OUTPUT)
	@touch $(%%_CHIMERA_OUTPUT)
