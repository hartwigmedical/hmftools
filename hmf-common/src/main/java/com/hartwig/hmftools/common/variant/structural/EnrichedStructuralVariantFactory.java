package com.hartwig.hmftools.common.variant.structural;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidyFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EnrichedStructuralVariantFactory {

    @NotNull
    public static List<EnrichedStructuralVariant> enrich(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers, @NotNull final List<StructuralVariant> variants) {

        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory =
                new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);

        final List<EnrichedStructuralVariant> result = Lists.newArrayList();
        for (final StructuralVariant variant : variants) {

            final ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);
            final ImmutableEnrichedStructuralVariantLeg.Builder startBuilder =
                    ImmutableEnrichedStructuralVariantLeg.builder().from(variant.start());
            final ImmutableEnrichedStructuralVariantLeg.Builder endBuilder =
                    ImmutableEnrichedStructuralVariantLeg.builder().from(variant.end());

            final List<StructuralVariantLegPloidy> ploidies = ploidyFactory.create(variant, copyNumbers);
            if (!ploidies.isEmpty()) {
                double roundedPloidy = round(ploidies.get(0).averageImpliedPloidy());
                builder.ploidy(roundedPloidy);
            }

            if (ploidies.size() == 2) {
                final StructuralVariantLegPloidy start = ploidies.get(0);
                final StructuralVariantLegPloidy end = ploidies.get(1);

                startBuilder.adjustedAlleleFrequency(round(adjustedVAF(purityAdjuster, start)));
                startBuilder.adjustedCopyNumber(round(adjustedCopyNumber(start)));
                startBuilder.adjustedCopyNumberChange(round(adjustedCopyNumberChange(start)));

                endBuilder.adjustedAlleleFrequency(round(adjustedVAF(purityAdjuster, end)));
                endBuilder.adjustedCopyNumber(round(adjustedCopyNumber(end)));
                endBuilder.adjustedCopyNumberChange(round(adjustedCopyNumberChange(end)));
            }

            result.add(builder.start(startBuilder.build()).end(endBuilder.build()).build());

        }

        return result;
    }

    @Nullable
    private static Double round(@Nullable Double value) {
        return value == null ? null : Math.round(value * 1000d) / 1000d;
    }

    @Nullable
    private static Double adjustedVAF(@NotNull final PurityAdjuster purityAdjuster, @NotNull final StructuralVariantLegPloidy ploidy) {
        final Double adjustedCopyNumber = adjustedCopyNumber(ploidy);
        return adjustedCopyNumber == null ? null : purityAdjuster.purityAdjustedVAF(ploidy.chromosome(), adjustedCopyNumber, ploidy.vaf());
    }

    @Nullable
    private static Double adjustedCopyNumber(@NotNull final StructuralVariantLegPloidy ploidy) {

        if (ploidy.orientation() == 1) {
            return ploidy.leftCopyNumber().orElse(null);
        } else {
            return ploidy.rightCopyNumber().orElse(null);
        }
    }

    @Nullable
    private static Double adjustedCopyNumberChange(@NotNull final StructuralVariantLegPloidy ploidy) {

        if (ploidy.leftCopyNumber().isPresent() && ploidy.rightCopyNumber().isPresent()) {
            return ploidy.orientation() == 1
                    ? ploidy.leftCopyNumber().get() - ploidy.rightCopyNumber().get()
                    : ploidy.rightCopyNumber().get() - ploidy.leftCopyNumber().get();
        }

        return null;
    }

}
