package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidyFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EnrichedStructuralVariantFactory {

    @NotNull
    public static List<EnrichedStructuralVariant> enrich(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers, @NotNull final List<StructuralVariant> variants) {

        final StructuralVariantPloidyFactory<PurpleCopyNumber> ploidyFactory =
                new StructuralVariantPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);

        final List<EnrichedStructuralVariant> result = Lists.newArrayList();
        for (final StructuralVariant variant : variants) {

            final ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);

            final List<StructuralVariantPloidy> ploidies = ploidyFactory.create(Collections.singletonList(variant), copyNumbers);
            if (!ploidies.isEmpty()) {
                double roundedPloidy = round(ploidies.get(0).averageImpliedPloidy());
                builder.ploidy(roundedPloidy);
            }

            if (ploidies.size() == 2) {
                final StructuralVariantPloidy start = ploidies.get(0);
                final StructuralVariantPloidy end = ploidies.get(1);

                builder.adjustedStartAF(round(adjustedVAF(purityAdjuster, start)));
                builder.adjustedStartCopyNumber(round(adjustedCopyNumber(start)));
                builder.adjustedStartCopyNumberChange(round(adjustedCopyNumberChange(start)));
                builder.adjustedEndAF(round(adjustedVAF(purityAdjuster, end)));
                builder.adjustedEndCopyNumber(round(adjustedCopyNumber(end)));
                builder.adjustedEndCopyNumberChange(round(adjustedCopyNumberChange(end)));
            }

            result.add(builder.build());

        }

        return result;
    }

    @Nullable
    private static Double round(@Nullable Double value) {
        return value == null ? null : Math.round(value * 1000d) / 1000d;
    }

    @Nullable
    private static Double adjustedVAF(@NotNull final PurityAdjuster purityAdjuster, @NotNull final StructuralVariantPloidy ploidy) {
        final Double adjustedCopyNumber = adjustedCopyNumber(ploidy);
        return adjustedCopyNumber == null ? null : purityAdjuster.purityAdjustedVAF(ploidy.chromosome(), adjustedCopyNumber, ploidy.vaf());
    }

    @Nullable
    private static Double adjustedCopyNumber(@NotNull final StructuralVariantPloidy ploidy) {

        if (ploidy.orientation() == 1) {
            return ploidy.leftCopyNumber().orElse(null);
        } else {
            return ploidy.rightCopyNumber().orElse(null);
        }
    }

    @Nullable
    private static Double adjustedCopyNumberChange(@NotNull final StructuralVariantPloidy ploidy) {

        if (ploidy.leftCopyNumber().isPresent() && ploidy.rightCopyNumber().isPresent()) {
            return ploidy.orientation() == 1
                    ? ploidy.leftCopyNumber().get() - ploidy.rightCopyNumber().get()
                    : ploidy.rightCopyNumber().get() - ploidy.leftCopyNumber().get();
        }

        return null;
    }

}
