package com.hartwig.hmftools.patientreporter.variants;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.and;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.jetbrains.annotations.NotNull;

class ConsequenceRule {

    private static final List<VariantConsequence> ACTIONABLE_CONSEQUENCES = Lists.newArrayList(
            VariantConsequence.TRANSCRIPT_ABLATION, VariantConsequence.TRANSCRIPT_AMPLIFICATION,
            VariantConsequence.SPLICE_ACCEPTOR_VARIANT, VariantConsequence.SPLICE_DONOR_VARIANT,
            VariantConsequence.SPLICE_REGION_VARIANT, VariantConsequence.STOP_GAINED, VariantConsequence.STOP_LOST,
            VariantConsequence.INCOMPLETE_TERMINAL_CODING_VARIANT, VariantConsequence.INITIATOR_CODON_VARIANT,
            VariantConsequence.START_LOST, VariantConsequence.FRAMESHIFT_VARIANT, VariantConsequence.INFRAME_INSERTION,
            VariantConsequence.INFRAME_DELETION, VariantConsequence.MISSENSE_VARIANT);

    @NotNull
    private final Slicer hmfSlicingRegion;

    ConsequenceRule(@NotNull final Slicer hmfSlicingRegion) {
        this.hmfSlicingRegion = hmfSlicingRegion;
    }

    @NotNull
    List<SomaticVariant> apply(@NotNull List<SomaticVariant> variants) {
        Predicate<SomaticVariant> consequenceRule = and(isIncludedIn(hmfSlicingRegion), hasActionableConsequence());

        return filter(variants, consequenceRule);
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull Slicer slicer) {
        return slicer::includes;
    }

    @NotNull
    private static Predicate<SomaticVariant> hasActionableConsequence() {
        return variant -> {
            for (VariantConsequence consequence : ACTIONABLE_CONSEQUENCES) {
                if (variant.hasConsequence(consequence)) {
                    return true;
                }
            }
            return false;
        };
    }
}
