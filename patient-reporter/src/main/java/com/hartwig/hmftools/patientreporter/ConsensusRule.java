package com.hartwig.hmftools.patientreporter;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.and;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.inCOSMIC;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.inDBSNP;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.not;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.or;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.withMinCallers;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.withType;

import java.util.List;
import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;

import org.jetbrains.annotations.NotNull;

class ConsensusRule {

    @NotNull
    private final Slicer highConfidenceRegion;
    @NotNull
    private final Slicer cpctSlicingRegion;

    ConsensusRule(@NotNull final Slicer highConfidenceRegion, @NotNull final Slicer cpctSlicingRegion) {
        this.highConfidenceRegion = highConfidenceRegion;
        this.cpctSlicingRegion = cpctSlicingRegion;
    }

    @NotNull
    List<SomaticVariant> apply(@NotNull List<SomaticVariant> variants) {
        Predicate<SomaticVariant> SNPRule = and(withType(VariantType.SNP),
                or(withMinCallers(3), isIncludedIn(cpctSlicingRegion),
                        and(withMinCallers(2), isIncludedIn(highConfidenceRegion), inCOSMIC(), not(inDBSNP()))));
        Predicate<SomaticVariant> IndelRule = and(withType(VariantType.INDEL),
                or(withMinCallers(2), isIncludedIn(cpctSlicingRegion)));

        return filter(variants, or(SNPRule, IndelRule));
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull Slicer slicer) {
        return slicer::includes;
    }
}
