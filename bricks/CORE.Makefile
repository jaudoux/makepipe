# GENERAL VARIABLES
%%_MAPPING_DIRECTORY = mapping
%%_CHIMERAS_DIRECTORY = chimeras

# GENERAL RULES
$(%%_MAPPING_DIRECTORY):
	mkdir -p $(%%_MAPPING_DIRECTORY)

$(%%_CHIMERAS_DIRECTORY):
	mkdir -p $(%%_CHIMERAS_DIRECTORY)

# DEFINING NEW IMPLICIT RULES
%.bam: %.sam
	samtools view -Sub $< | samtools sort -@ 10 -m 2G - $*
	samtools index $@
