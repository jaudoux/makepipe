## THIS FILE HAS BEEN GENERATED BY 'makePipe v0.0' from pipeline.yml

all: $(lastword $(MAKEFILE_LIST))
	$(MAKE) -j 4 counts chimCT_test1 chimCT_test2 cractools_test1 cractools_test2 crac_test1 crac_test2
.PHONY: all

$(lastword $(MAKEFILE_LIST)): pipeline.yml
	./makepipe build $< > $@

##############################
#          VARIABLES         #
##############################

# Process: counts, constructed from BRICK FEATURECOUNTS
counts_ANNOTATION_FILE = /data/annotations/human/Homo_sapiens.GRCh38.77.gtf
counts_INPUT_ALIGNMENTS = mapping/crac/test1.bam mapping/crac/test2.bam
counts_OUTPUT_FILE = quantification/gene_counts.tsv
counts_BINARY = featureCounts # Binary used to run FEATURECOUNTS
counts_OPTIONS = # Options that can be specified to FEATURECOUNTS (run featurecounts -h for more informations


# Process: chimCT_test1, constructed from BRICK CHIMCT
chimCT_ANNOTATION_FILE = /data/annotations/human/ensembl.gff
chimCT_BINARY = chimCT # Binary used to call CHIMCT
chimCT_OPTIONS = # Options provided to CHIMCT

chimCT_test1_SAM_FILE = mapping/crac/test1.bam
chimCT_test1_OUTPUT_FILE = chimeras/chimCT/test1-chimeras.csv

# Process: chimCT_test2, constructed from BRICK CHIMCT

chimCT_test2_SAM_FILE = mapping/crac/test2.bam
chimCT_test2_OUTPUT_FILE = chimeras/chimCT/test2-chimeras.csv

# Process: cractools_test1, constructed from BRICK CRACTOOLS_EXTRACT
cractools_OPTIONS = --coverless-splices # Options passed to CRACTOOLS_EXTRACT
cractools_BINARY = cractools-extract # Binary used to call CRACTOOLS_EXTRACT

cractools_test1_SAM_FILE = mapping/crac/test1.bam
cractools_test1_CHIMERA_OUTPUT = chimeras/cractools/test1-chimeras.tsv
cractools_test1_SPLICE_OUTPUT = splices/cractools/test1-splice.bed
cractools_test1_MUTATION_OUTPUT = mutations/cractools/test1.vcf

# Process: cractools_test2, constructed from BRICK CRACTOOLS_EXTRACT

cractools_test2_SAM_FILE = mapping/crac/test2.bam
cractools_test2_CHIMERA_OUTPUT = chimeras/cractools/test2-chimeras.tsv
cractools_test2_SPLICE_OUTPUT = splices/cractools/test2-splice.bed
cractools_test2_MUTATION_OUTPUT = mutations/cractools/test2.vcf

# Process: crac_test1, constructed from BRICK CRAC
crac_INDEX = /data/indexes/crac/GRCh38
crac_BINARY = /data/projects/crac-dev/src/crac
crac_KMER_LENGTH = 22
crac_OPTIONS = --no-ambiguity
crac_VERSION = $(shell $(crac_BINARY) -v | head -1 | cut -d " " -f 3)
crac_OPTIONS += -i $(crac_INDEX) -k $(crac_KMER_LENGTH) --bam
crac_OUTPUT_EXTENSION = $(suffix $(crac_test1_OUTPUT_FILE))

crac_test1_READS = raw/test1_1.fastq raw/test1_2.fastq
crac_test1_OUTPUT_FILE = mapping/crac/test1.bam

# Process: crac_test2, constructed from BRICK CRAC

crac_test2_READS = raw/test2_1.fastq raw/test2_2.fastq
crac_test2_OUTPUT_FILE = mapping/crac/test2.bam

##############################
#            RULES           #
##############################

# Process: counts, constructed from BRICK FEATURECOUNTS
counts: $(counts_OUTPUT_FILE)
counts_clean:
	-rm $(basename $(counts_OUTPUT_FILE))*
	-rmdir -p --ignore-fail-on-non-empty $(dir $(counts_OUTPUT_FILE))
$(counts_OUTPUT_FILE): $(counts_INPUT_ALIGNMENTS) $(counts_ANNOTATION_FILE)
	@mkdir -p $(dir $@)
	$(counts_BINARY) $(counts_OPTIONS) -a $(counts_ANNOTATION_FILE) -o $@ $(counts_INPUT_ALIGNMENTS) > $@

# Process: chimCT_test1, constructed from BRICK CHIMCT
.PHONY: chimCT_test1
chimCT_test1: $(chimCT_test1_OUTPUT_FILE) 
chimCT_test1_clean:
	-rm $(basename $(chimCT_test1_OUTPUT_FILE))*
	-rmdir -p --ignore-fail-on-non-empty $(dir $(chimCT_test1_OUTPUT_FILE))
$(chimCT_test1_OUTPUT_FILE): $(chimCT_test1_SAM_FILE)
	@mkdir -p $(dir $@)
	$(chimCT_BINARY) -n $(basename $(notdir $<)) -g $(chimCT_ANNOTATION_FILE) -s $< $(chimCT_OPTIONS) > $@ 2> $(@:$(suffix $@)=.log)

# Process: chimCT_test2, constructed from BRICK CHIMCT
.PHONY: chimCT_test2
chimCT_test2: $(chimCT_test2_OUTPUT_FILE) 
chimCT_test2_clean:
	-rm $(basename $(chimCT_test2_OUTPUT_FILE))*
	-rmdir -p --ignore-fail-on-non-empty $(dir $(chimCT_test2_OUTPUT_FILE))
$(chimCT_test2_OUTPUT_FILE): $(chimCT_test2_SAM_FILE)
	@mkdir -p $(dir $@)
	$(chimCT_BINARY) -n $(basename $(notdir $<)) -g $(chimCT_ANNOTATION_FILE) -s $< $(chimCT_OPTIONS) > $@ 2> $(@:$(suffix $@)=.log)

# Process: cractools_test1, constructed from BRICK CRACTOOLS_EXTRACT
cractools_test1: $(cractools_test1_SPLICE_OUTPUT)
cractools_test1_clean:
	-rm $(cractools_test1_MUTATION_OUTPUT)
	-rm $(cractools_test1_SPLICE_OUTPUT)
	-rm $(cractools_test1_CHIMERA_OUTPUT)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(cractools_test1_MUTATION_OUTPUT))
	-rmdir -p --ignore-fail-on-non-empty $(dir $(cractools_test1_SPLICE_OUTPUT))
	-rmdir -p --ignore-fail-on-non-empty $(dir $(cractools_test1_CHIMERA_OUTPUT))
$(cractools_test1_MUTATION_OUTPUT): $(cractools_test1_SPLICE_OUTPUT)
$(cractools_test1_CHIMERA_OUTPUT): $(cractools_test1_SPLICE_OUTPUT)
$(cractools_test1_SPLICE_OUTPUT): $(cractools_test1_SAM_FILE)
	@mkdir -p $(dir $(cractools_test1_CHIMERA_OUTPUT) $(cractools_test1_SPLICE_OUTPUT) $(cractools_test1_MUTATION_OUTPUT))
	$(cractools_BINARY) $(cractools_OPTIONS) $< -s $@ -m $(cractools_test1_MUTATION_OUTPUT) -c $(cractools_test1_CHIMERA_OUTPUT)
	@touch $(cractools_test1_MUTATION_OUTPUT)
	@touch $(cractools_test1_CHIMERA_OUTPUT)

# Process: cractools_test2, constructed from BRICK CRACTOOLS_EXTRACT
cractools_test2: $(cractools_test2_SPLICE_OUTPUT)
cractools_test2_clean:
	-rm $(cractools_test2_MUTATION_OUTPUT)
	-rm $(cractools_test2_SPLICE_OUTPUT)
	-rm $(cractools_test2_CHIMERA_OUTPUT)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(cractools_test2_MUTATION_OUTPUT))
	-rmdir -p --ignore-fail-on-non-empty $(dir $(cractools_test2_SPLICE_OUTPUT))
	-rmdir -p --ignore-fail-on-non-empty $(dir $(cractools_test2_CHIMERA_OUTPUT))
$(cractools_test2_MUTATION_OUTPUT): $(cractools_test2_SPLICE_OUTPUT)
$(cractools_test2_CHIMERA_OUTPUT): $(cractools_test2_SPLICE_OUTPUT)
$(cractools_test2_SPLICE_OUTPUT): $(cractools_test2_SAM_FILE)
	@mkdir -p $(dir $(cractools_test2_CHIMERA_OUTPUT) $(cractools_test2_SPLICE_OUTPUT) $(cractools_test2_MUTATION_OUTPUT))
	$(cractools_BINARY) $(cractools_OPTIONS) $< -s $@ -m $(cractools_test2_MUTATION_OUTPUT) -c $(cractools_test2_CHIMERA_OUTPUT)
	@touch $(cractools_test2_MUTATION_OUTPUT)
	@touch $(cractools_test2_CHIMERA_OUTPUT)

# Process: crac_test1, constructed from BRICK CRAC
.PHONY: crac_test1
crac_test1: $(crac_test1_OUTPUT_FILE)
crac_test1_clean:
	-rm $(basename $(crac_test1_OUTPUT_FILE))*
	-rmdir -p --ignore-fail-on-non-empty $(dir $(crac_test1_OUTPUT_FILE))
.SERIAL: $(crac_test1_OUTPUT_FILE)
$(crac_test1_OUTPUT_FILE): $(crac_test1_READS)
	mkdir -p $(dir $@)
	$(crac_BINARY) $(crac_OPTIONS) -r $(wordlist 1, 2, $^) -o - --summary $(@:$(crac_OUTPUT_EXTENSION)=-summary.txt) 2> $(@:$(crac_OUTPUT_EXTENSION)=.log) | samtools sort - $(@:%$(crac_OUTPUT_EXTENSION)=%)
	samtools index $@

# Process: crac_test2, constructed from BRICK CRAC
.PHONY: crac_test2
crac_test2: $(crac_test2_OUTPUT_FILE)
crac_test2_clean:
	-rm $(basename $(crac_test2_OUTPUT_FILE))*
	-rmdir -p --ignore-fail-on-non-empty $(dir $(crac_test2_OUTPUT_FILE))
.SERIAL: $(crac_test2_OUTPUT_FILE)
$(crac_test2_OUTPUT_FILE): $(crac_test2_READS)
	mkdir -p $(dir $@)
	$(crac_BINARY) $(crac_OPTIONS) -r $(wordlist 1, 2, $^) -o - --summary $(@:$(crac_OUTPUT_EXTENSION)=-summary.txt) 2> $(@:$(crac_OUTPUT_EXTENSION)=.log) | samtools sort - $(@:%$(crac_OUTPUT_EXTENSION)=%)
	samtools index $@

clean: counts_clean chimCT_test1_clean chimCT_test2_clean cractools_test1_clean cractools_test2_clean crac_test1_clean crac_test2_clean