package com.hartwig.hmftools.common.variant.confidenceregions;

import static com.hartwig.hmftools.common.variant.predicate.VariantFilter.filter;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.and;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.inCOSMIC;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.inDBSNP;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.not;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.or;
import static com.hartwig.hmftools.common.variant.predicate.VariantPredicates.withType;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;

import org.jetbrains.annotations.NotNull;

public class ConfidenceRegionsRule {

    @VisibleForTesting
    public static final String FILTER = "CONFIDENCE_REGIONS";

    @NotNull
    private final Predicate<SomaticVariant> confidenceRegionsFilter;

    @NotNull
    public static ConfidenceRegionsRule fromSlicers(@NotNull final Slicer highConfidenceRegion,
            @NotNull final Slicer extremeConfidenceRegion) {
        final Predicate<SomaticVariant> snpRule = and(withType(VariantType.SNP),
                or(isIncludedIn(extremeConfidenceRegion), and(isIncludedIn(highConfidenceRegion), or(inCOSMIC(), not(inDBSNP())))));

        final Predicate<SomaticVariant> indelRule = and(withType(VariantType.INDEL), isIncludedIn(extremeConfidenceRegion));

        return new ConfidenceRegionsRule(or(snpRule, indelRule));
    }

    @VisibleForTesting
    ConfidenceRegionsRule(@NotNull final Predicate<SomaticVariant> confidenceRegionsFilter) {
        this.confidenceRegionsFilter = confidenceRegionsFilter;
    }

    @NotNull
    public List<SomaticVariant> removeUnreliableVariants(@NotNull final List<SomaticVariant> variants) {
        return filter(variants, confidenceRegionsFilter);
    }

    @NotNull
    public List<SomaticVariant> updateFilterFlagForUnreliableVariants(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> updatedVariants = Lists.newArrayList();
        for (final SomaticVariant variant : variants) {
            final SomaticVariant.Builder newVariantBuilder = SomaticVariant.Builder.fromVariant(variant);
            if (!confidenceRegionsFilter.test(variant)) {
                if (VariantFilter.isPass(variant)) {
                    newVariantBuilder.filter(FILTER);
                } else {
                    newVariantBuilder.filter(variant.filter() + ";" + FILTER);
                }
            }
            updatedVariants.add(newVariantBuilder.build());
        }
        return updatedVariants;
    }

    @NotNull
    private static Predicate<SomaticVariant> isIncludedIn(@NotNull Slicer slicer) {
        return slicer::includes;
    }
}
