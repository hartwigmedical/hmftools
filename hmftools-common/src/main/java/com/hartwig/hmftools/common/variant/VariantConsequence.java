package com.hartwig.hmftools.common.variant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public enum VariantConsequence {
    TRANSCRIPT_ABLATION("transcript_ablation"),
    TRANSCRIPT_AMPLIFICATION("transcript_amplification"),
    SPLICE_ACCEPTOR_VARIANT("splice_acceptor_variant"),
    SPLICE_DONOR_VARIANT("splice_donor_variant"),
    SPLICE_REGION_VARIANT("splice_region_variant"),
    STOP_GAINED("stop_gained"),
    STOP_LOST("stop_lost"),
    INCOMPLETE_TERMINAL_CODING_VARIANT("incomplete_terminal_coding_variant"),
    INITIATOR_CODON_VARIANT("initiator_codon_variant"),
    START_LOST("start_lost"),
    FRAMESHIFT_VARIANT("frameshift_variant"),
    INFRAME_INSERTION("inframe_insertion"),
    INFRAME_DELETION("inframe_deletion"),
    MISSENSE_VARIANT("missense_variant"),
    OTHER(Strings.EMPTY);

    @NotNull
    private final String sequenceOntologyTerm;

    VariantConsequence(@NotNull final String sequenceOntologyTerm) {
        this.sequenceOntologyTerm = sequenceOntologyTerm;
    }

    @NotNull
    public String sequenceOntologyTerm() {
        return sequenceOntologyTerm;
    }
}
