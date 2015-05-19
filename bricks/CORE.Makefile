## GENERAL VARIABLES
%%_MAPPING_DIRECTORY 				= mapping
%%_CHIMERAS_DIRECTORY 			= chimeras
%%_SPLICES_DIRECTORY 				= splices
%%_MUTATIONS_DIRECTORY 			= mutations
%%_QUANTIFICATION_DIRECTORY = quantification

%%_DIRECTORIES = $(%%_MAPPING_DIRECTORY) $(%%_CHIMERAS_DIRECTORY) $(%%_SPLICES_DIRECTORY) $(%%_MUTATIONS_DIRECTORY) $(%%_QUANTIFICATION_DIRECTORY)

## MAKING DIRECTORY RULES
$(%%_DIRECTORIES):
	mkdir -p $@

#DEFINING NEW IMPLICIT RULES
%.bam: %.sam
	samtools view -Sub $< | samtools sort -@ 10 -m 2G - $*
	samtools index $@
