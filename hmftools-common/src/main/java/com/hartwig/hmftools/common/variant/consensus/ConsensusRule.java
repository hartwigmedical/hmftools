package com.hartwig.hmftools.common.variant.consensus;

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

import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class ConsensusRule {

    @NotNull
    private final Slicer highConfidenceRegion;
    @NotNull
    private final Slicer extremeConfidenceRegion;

    public ConsensusRule(@NotNull final Slicer highConfidenceRegion, @NotNull final Slicer extremeConfidenceRegion) {
        this.highConfidenceRegion = highConfidenceRegion;
        this.extremeConfidenceRegion = extremeConfidenceRegion;
    }

    @NotNull
    public List<SomaticVariant> apply(@NotNull List<SomaticVariant> variants) {
        final Predicate<SomaticVariant> snpRule = and(withType(VariantType.SNP),
                or(withMinCallers(3), isIncludedIn(extremeConfidenceRegion),
                        and(withMinCallers(2), isIncludedIn(highConfidenceRegion),
                                or(inCOSMIC(), not(inDBSNP())))));
        final Predicate<SomaticVariant> indelRule = and(withType(VariantType.INDEL),
                or(withMinCallers(2), isIncludedIn(extremeConfidenceRegion)));

        return filter(variants, or(snpRule, indelRule));
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull Slicer slicer) {
        return slicer::includes;
    }
}
