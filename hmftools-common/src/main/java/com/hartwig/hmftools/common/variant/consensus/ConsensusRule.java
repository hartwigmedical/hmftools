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
    private final Slicer giabHighConfidenceRegion;
    @NotNull
    private final Slicer cpctSlicingRegion;

    public ConsensusRule(@NotNull final Slicer giabHighConfidenceRegion, @NotNull final Slicer cpctSlicingRegion) {
        this.giabHighConfidenceRegion = giabHighConfidenceRegion;
        this.cpctSlicingRegion = cpctSlicingRegion;
    }

    @NotNull
    public List<SomaticVariant> apply(@NotNull List<SomaticVariant> variants) {
        final Predicate<SomaticVariant> snpRule = and(withType(VariantType.SNP),
                or(withMinCallers(3), isIncludedIn(cpctSlicingRegion),
                        and(withMinCallers(2), isIncludedIn(giabHighConfidenceRegion),
                                or(inCOSMIC(), not(inDBSNP())))));
        final Predicate<SomaticVariant> indelRule = and(withType(VariantType.INDEL),
                or(withMinCallers(2), isIncludedIn(cpctSlicingRegion)));

        return filter(variants, or(snpRule, indelRule));
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull Slicer slicer) {
        return slicer::includes;
    }
}
