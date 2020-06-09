package com.hartwig.hmftools.common.variant.structural;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegCopyNumberChangeFactory;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegCopyNumberFactory;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidyFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CopyNumberEnrichedStructuralVariantFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;

    public CopyNumberEnrichedStructuralVariantFactory(@NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Multimap<Chromosome, PurpleCopyNumber> copyNumbers) {
        this.purityAdjuster = purityAdjuster;
        this.copyNumbers = copyNumbers;
    }

    @NotNull
    public List<EnrichedStructuralVariant> enrich(@NotNull List<StructuralVariant> variants) {
        final StructuralVariantLegCopyNumberChangeFactory changeFactory =
                new StructuralVariantLegCopyNumberChangeFactory(purityAdjuster, copyNumbers, variants);

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

            @Nullable
            final StructuralVariantLeg endLeg = variant.end();
            ImmutableEnrichedStructuralVariantLeg.Builder endBuilder = createBuilder(endLeg);

            List<StructuralVariantLegPloidy> ploidies = ploidyFactory.create(variant, copyNumbers);
            if (!ploidies.isEmpty()) {
                // The implied ploidy should be equal between start and end, so doesn't matter what we pick.
                builder.junctionCopyNumber((ploidies.get(0).averageImpliedPloidy()));

                StructuralVariantLegPloidy startPloidy = ploidies.get(0);
                StructuralVariantLegPloidy endPloidy = ploidies.size() <= 1 ? null : ploidies.get(1);

                startBuilder.adjustedAlleleFrequency((startPloidy.adjustedVaf()));
                startBuilder.adjustedCopyNumber((startPloidy.adjustedCopyNumber()));
                startBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(startPloidy)));

                if (endPloidy != null) {
                    assert endBuilder != null;
                    endBuilder.adjustedAlleleFrequency((endPloidy.adjustedVaf()));
                    endBuilder.adjustedCopyNumber((endPloidy.adjustedCopyNumber()));
                    endBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(endPloidy)));
                }
            } else {
                // Can't always get ploidies (if no vaf for example) but we can still get copy number info
                final StructuralVariantLegCopyNumber startCopyNumber = copyNumberFactory.create(variant.start(), copyNumbers);
                startBuilder.adjustedCopyNumber((startCopyNumber.adjustedCopyNumber()));
                startBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(startCopyNumber)));

                // Lacking anything else, inferred variants can use copy number change as ploidy
                if (variant.type() == StructuralVariantType.INF) {
                    builder.junctionCopyNumber((changeFactory.copyNumberChange(startCopyNumber)));
                }

                if (endLeg != null) {
                    assert endBuilder != null;
                    final StructuralVariantLegCopyNumber endCopyNumber = copyNumberFactory.create(endLeg, copyNumbers);
                    endBuilder.adjustedCopyNumber((endCopyNumber.adjustedCopyNumber()));
                    endBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(endCopyNumber)));
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
}
