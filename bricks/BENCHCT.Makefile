## GENERAL
%%_OUTPUT_FILE=undef
%%_CONFIG_FILE=benchCT.yaml
%%_CONFIG_FILE_DEPENDENCIES=
%%_BINARY=benchCT
%%_NB_THREADS=1
%%_DEPENDENCIES=

## TRUTH FILES
%%_CHECK_INFO_FILE=
%%_CHECK_SPLICE_FILE=
%%_CHECK_CHIMERA_FILE= 
%%_CHECK_MUTATION_FILE=
%%_CHECK_TRANSCRIPT_FILE=

%%_DEPENDENCIES += $(%%_CHECK_MAPPING_FILE) $(%%_CHECK_CHIMERAS_FILE) $(%%_CHECK_SPLICES_FILE) $(%%_CHECK_MUTATIONS_FILE) $(%%_CHECK_TRANSCRIPT_FILE) $(%%_CHECK_ERROR_FILE)

%%_SOFTWARES_CONFIG=

.PHONY: %%
%%: $(%%_OUTPUT_FILE)

%%_clean:
	-rm $(%%_OUTPUT_FILE)
	-rmdir -p --ignore-fail-on-non-empty $(dir $(%%_OUTPUT_FILE))

$(%%_OUTPUT_FILE): $(%%_CONFIG_FILE)
	$(%%_BINARY) -v $< > $@ -p $(%%_NB_THREADS)

$(%%_CONFIG_FILE): $(%%_DEPENDENCIES)
	@echo "---" > $@
	@echo "checker:" >> $@
	@echo "    files:" >> $@
ifneq ($(%%_CHECK_INFO_FILE),)
	@echo "        infos:     $(%%_CHECK_INFO_FILE)" >> $@
endif
ifneq ($(%%_CHECK_SPLICE_FILE),)
	@echo "        splices:     $(%%_CHECK_SPLICE_FILE)" >> $@
endif
ifneq ($(%%_CHECK_MUTATION_FILE),)
	@echo "        mutations:   $(%%_CHECK_MUTATION_FILE)" >> $@
endif
ifneq ($(%%_CHECK_CHIMERA_FILE),)
	@echo "        chimeras:    $(%%_CHECK_CHIMERA_FILE)" >> $@
endif
ifneq ($(%%_CHECK_TRANSCRIPT_FILE),)
	@echo "        annotations:    $(%%_CHECK_TRANSCRIPT_FILE)" >> $@
endif
	@echo "output:" >> $@
	@echo "    statistics:" >> $@
	@echo "        - sensitivity" >> $@
	@echo "        - accuracy" >> $@
	@echo "        - true-positives" >> $@
	@echo "        - false-positives" >> $@
	@echo "        - nb-elements" >> $@
	@echo "    nb_decimals: 4" >> $@
	echo -E "softwares: " '$(%%_SOFTWARES_CONFIG)' >> $@
