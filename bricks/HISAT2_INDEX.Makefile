## INPUT VARS:
%%_FASTA_FILES = undef

## OPTIONS
%%_INDEX_PREFIX = undef
%%_BINARY = hisat2-build
%%_OPTIONS = #

%%_COMMA = ,
%%_SPACE = $(empty) $(empty)

%%: $(%%_INDEX_PREFIX)

.PHONY: $(%%_INDEX_PREFIX)
$(%%_INDEX_PREFIX): $(%%_INDEX_PREFIX).1.ht2

%%_clean:
	-rm -rf $(%%_INDEX_PREFIX).*

$(%%_INDEX_PREFIX).1.ht2:	$(%%_FASTA_FILES)
	mkdir -p $(dir $@)
	$(%%_BINARY) $(%%_OPTIONS) $(subst $(%%_SPACE),$(%%_COMMA),$^) $(%%_INDEX_PREFIX)
