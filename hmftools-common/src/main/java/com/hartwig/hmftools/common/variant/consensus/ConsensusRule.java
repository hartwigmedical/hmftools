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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;

import org.jetbrains.annotations.NotNull;

public class ConsensusRule {

    @VisibleForTesting
    static final String CONSENSUS_FILTERED = "CALLER_CONSENSUS";

    @NotNull
    private final Predicate<SomaticVariant> consensusRuleFilter;

    @NotNull
    public static ConsensusRule fromSlicers(@NotNull final Slicer highConfidenceRegion,
            @NotNull final Slicer extremeConfidenceRegion) {
        final Predicate<SomaticVariant> snpRule = and(withType(VariantType.SNP),
                or(withMinCallers(3), isIncludedIn(extremeConfidenceRegion),
                        and(withMinCallers(2), isIncludedIn(highConfidenceRegion), or(inCOSMIC(), not(inDBSNP())))));
        final Predicate<SomaticVariant> indelRule = and(withType(VariantType.INDEL),
                or(withMinCallers(2), isIncludedIn(extremeConfidenceRegion)));

        return new ConsensusRule(or(snpRule, indelRule));
    }

    @VisibleForTesting
    ConsensusRule(@NotNull final Predicate<SomaticVariant> consensusRuleFilter) {
        this.consensusRuleFilter = consensusRuleFilter;
    }

    @NotNull
    public List<SomaticVariant> removeUnreliableVariants(@NotNull final List<SomaticVariant> variants) {
        return filter(variants, consensusRuleFilter);
    }

    @NotNull
    public List<SomaticVariant> updateFilterFlagForUnreliableVariants(@NotNull final List<SomaticVariant> variants) {
        final List<SomaticVariant> updatedVariants = Lists.newArrayList();
        for (final SomaticVariant variant : variants) {
            final SomaticVariant.Builder newVariantBuilder = SomaticVariant.Builder.fromVariant(variant);
            if (!consensusRuleFilter.test(variant)) {
                if (VariantFilter.isPass(variant)) {
                    newVariantBuilder.filter(CONSENSUS_FILTERED);
                } else {
                    newVariantBuilder.filter(variant.filter() + ";" + CONSENSUS_FILTERED);
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
