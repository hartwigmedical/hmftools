package com.hartwig.hmftools.common.variant.structural;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidyFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberEnrichedStructuralVariantFactory {

    private final PurityAdjuster purityAdjuster;
    private final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;

    public CopyNumberEnrichedStructuralVariantFactory(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Multimap<Chromosome, PurpleCopyNumber> copyNumbers) {
        this.purityAdjuster = purityAdjuster;
        this.copyNumbers = copyNumbers;
    }

    @NotNull
    public List<EnrichedStructuralVariant> enrich(@NotNull final List<StructuralVariant> variants) {
        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory =
                new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
        final StructuralVariantLegCopyNumberFactory<PurpleCopyNumber> copyNumberFactory =
                new StructuralVariantLegCopyNumberFactory<>(PurpleCopyNumber::averageTumorCopyNumber);

        final List<EnrichedStructuralVariant> result = Lists.newArrayList();
        for (final StructuralVariant variant : variants) {


            ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);
            ImmutableEnrichedStructuralVariantLeg.Builder startBuilder = createBuilder(variant.start());
            // Every SV should have a start, while end is optional.
            assert startBuilder != null;

            @Nullable final StructuralVariantLeg endLeg = variant.end();
            ImmutableEnrichedStructuralVariantLeg.Builder endBuilder = createBuilder(endLeg);

            List<StructuralVariantLegPloidy> ploidies = ploidyFactory.create(variant, copyNumbers);
            if (!ploidies.isEmpty()) {
                // The implied ploidy should be equal between start and end, so doesn't matter what we pick.
                builder.ploidy(round(ploidies.get(0).averageImpliedPloidy()));

                StructuralVariantLegPloidy startPloidy = ploidies.get(0);
                StructuralVariantLegPloidy endPloidy = ploidies.size() <= 1 ? null : ploidies.get(1);

                startBuilder.adjustedAlleleFrequency(round(startPloidy.adjustedVaf()));
                startBuilder.adjustedCopyNumber(round(startPloidy.adjustedCopyNumber()));
                startBuilder.adjustedCopyNumberChange(round(startPloidy.adjustedCopyNumberChange()));

                if (endPloidy != null) {
                    assert endBuilder != null;
                    endBuilder.adjustedAlleleFrequency(round(endPloidy.adjustedVaf()));
                    endBuilder.adjustedCopyNumber(round(endPloidy.adjustedCopyNumber()));
                    endBuilder.adjustedCopyNumberChange(round(endPloidy.adjustedCopyNumberChange()));
                }
            } else {
                // Can't always get plodies (if no vaf for example) but we can still get copy number info
                final StructuralVariantLegCopyNumber startCopyNumber = copyNumberFactory.create(variant.start(), copyNumbers);
                startBuilder.adjustedCopyNumber(round(startCopyNumber.adjustedCopyNumber()));
                startBuilder.adjustedCopyNumberChange(round(startCopyNumber.adjustedCopyNumberChange()));

                if (endLeg != null) {
                    assert endBuilder != null;
                    final StructuralVariantLegCopyNumber endCopyNumber = copyNumberFactory.create(endLeg, copyNumbers);
                    endBuilder.adjustedCopyNumber(round(endCopyNumber.adjustedCopyNumber()));
                    endBuilder.adjustedCopyNumberChange(round(endCopyNumber.adjustedCopyNumberChange()));
                }
            }


            result.add(builder.start(startBuilder.build()).end(endBuilder == null ? null : endBuilder.build()).build());
        }

        return result;
    }

    @Nullable
    private ImmutableEnrichedStructuralVariantLeg.Builder createBuilder(@Nullable StructuralVariantLeg leg) {
        if (leg == null) {
            return null;
        }

        return ImmutableEnrichedStructuralVariantLeg.builder().from(leg);
    }

    private static double round(double value) {
        return Math.round(value * 1000d) / 1000d;
    }
}
