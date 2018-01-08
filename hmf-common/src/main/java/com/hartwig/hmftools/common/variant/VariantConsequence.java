package com.hartwig.hmftools.common.variant;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public enum VariantConsequence {
    // KODU: See also http://sequenceontology.org
    TRANSCRIPT(false, "transcript"),
    NON_CODING_EXON_VARIANT(false, "non_coding_exon_variant", "non_coding_transcript_exon_variant"),
    INTRON_VARIANT(false, "intron_variant"),
    INTRAGENIC_VARIANT(false, "intragenic_variant"),
    SEQUENCE_FEATURE(false, "sequence_feature"),
    SYNONYMOUS_VARIANT(false, "synonymous_variant", "stop_retained_variant"),
    UTR_VARIANT(false,
            "UTR_variant",
            "3_prime_UTR_variant",
            "5_prime_UTR_variant",
            "5_prime_UTR_premature_start_codon_gain_variant",
            "5_prime_UTR_truncation",
            "3_prime_UTR_truncation"),
    REGULATORY_REGION_VARIANT(false, "regulatory_region_variant", "TF_binding_site_variant"),
    INITIATOR_CODON_VARIANT(false, "initiator_codon_variant"),
    EXON_LOSS_VARIANT(false, "exon_loss_variant", "exon_loss"),
    NON_CANONICAL_START_CODON(false, "non_canonical_start_codon"),
    TRANSCRIPT_ABLATION(true, "transcript_ablation"),
    TRANSCRIPT_AMPLIFICATION(true, "transcript_amplification"),
    SPLICE_ACCEPTOR_VARIANT(true, "splice_acceptor_variant"),
    SPLICE_DONOR_VARIANT(true, "splice_donor_variant"),
    SPLICE_REGION_VARIANT(true, "splice_region_variant", "exonic_splice_region_variant", "non_coding_transcript_splice_region_variant"),
    STOP_GAINED(true, "stop_gained"),
    STOP_LOST(true, "stop_lost"),
    START_LOST(true, "start_lost"),
    FRAMESHIFT_VARIANT(true,
            "frameshift_variant",
            "frame_restoring_variant",
            "frameshift_elongation",
            "frameshift_truncation",
            "minus_1_frameshift_variant",
            "minus_2_frameshift_variant",
            "plus_1_frameshift_variant",
            "plus_2_frameshift_variant"),
    INFRAME_INSERTION(true, "inframe_insertion", "conservative_inframe_insertion", "disruptive_inframe_insertion"),
    INFRAME_DELETION(true, "inframe_deletion", "conservative_inframe_deletion", "disruptive_inframe_deletion"),
    MISSENSE_VARIANT(true,
            "missense_variant",
            "conservative_missense_variant",
            "non_conservative_missense_variant",
            "rare_amino_acid_variant",
            "pyrrolysine_loss",
            "selenocysteine_loss"),
    STRUCTURAL_INTERACTION_VARIANT(false, "structural_interaction_variant"),
    OTHER(false, Strings.EMPTY);

    @NotNull
    private final String parentSequenceOntologyTerm;
    @NotNull
    private final List<String> sequenceOntologySubTerms;
    private final boolean actionable;

    VariantConsequence(final boolean actionable, @NotNull final String parentSequenceOntologyTerm,
            @NotNull final String... sequenceOntologySubTerms) {
        this.parentSequenceOntologyTerm = parentSequenceOntologyTerm;
        this.sequenceOntologySubTerms = Lists.newArrayList(sequenceOntologySubTerms);
        this.actionable = actionable;
    }

    public boolean isParentTypeOf(@NotNull final String annotation) {
        return annotation.equals(parentSequenceOntologyTerm) || sequenceOntologySubTerms.contains(annotation);
    }

    @NotNull
    public String readableSequenceOntologyTerm() {
        return parentSequenceOntologyTerm.replace("_", " ");
    }

    public boolean isActionable() {
        return actionable;
    }

    @NotNull
    public static final List<VariantConsequence> ACTIONABLE_CONSEQUENCES =
            Stream.of(VariantConsequence.values()).filter(VariantConsequence::isActionable).collect(Collectors.toList());
}
