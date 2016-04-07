%%_BINARY=simCT
%%_OUTPUT_DIR=undef
%%_GENOME_DIR=undef
%%_ANNOTATIONS=undef
%%_OPTIONS=
%%_VCF_FILE=

%%_INFO_FILE=$(%%_OUTPUT_DIR)/info.txt
%%_READS_FILE=$(addprefix $(%%_OUTPUT_DIR)/reads_,$(addsuffix .fastq.gz, 1 2))

ifneq ($(%%_VCF_FILE),)
%%_OPTIONS+=--vcf-file $(%%_VCF_FILE)
endif

.PHONY: %%
%%: $(%%_OUTPUT_DIR)/info.txt

$(%%_READS_FILE): $(%%_INFO_FILE)

$(%%_INFO_FILE): $(%%_VCF_FILE)
	$(%%_BINARY) -g $(%%_GENOME_DIR) -o $(%%_OUTPUT_DIR) -a $(%%_ANNOTATIONS) $(%%_OPTIONS)
