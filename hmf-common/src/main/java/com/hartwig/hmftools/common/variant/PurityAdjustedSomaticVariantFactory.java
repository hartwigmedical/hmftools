package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.jetbrains.annotations.NotNull;

public class PurityAdjustedSomaticVariantFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;

    public PurityAdjustedSomaticVariantFactory(@NotNull PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        this.purityAdjuster = purityAdjuster;
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
    }

    public PurityAdjustedSomaticVariantFactory(@NotNull PurityAdjuster purityAdjuster,
            @NotNull final Multimap<String, PurpleCopyNumber> copyNumbers) {
        this.purityAdjuster = purityAdjuster;
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
    }

    @NotNull
    public List<PurityAdjustedSomaticVariant> create(@NotNull List<SomaticVariant> variants) {
        final List<PurityAdjustedSomaticVariant> result = Lists.newArrayList();

        for (SomaticVariant variant : variants) {
            final ImmutablePurityAdjustedSomaticVariantImpl.Builder builder = ImmutablePurityAdjustedSomaticVariantImpl.builder()
                    .from(variant)
                    .adjustedCopyNumber(0)
                    .adjustedVAF(0)
                    .ploidy(0)
                    .clonality(Clonality.UNKNOWN)
                    .lossOfHeterozygosity(false);

            copyNumberSelector.select(variant).ifPresent(x -> purityAdjustment(x, variant, builder));
            result.add(builder.build());
        }

        return result;
    }

    private void purityAdjustment(@NotNull final PurpleCopyNumber copyNumber, @NotNull final AllelicDepth depth,
            @NotNull final ImmutablePurityAdjustedSomaticVariantImpl.Builder builder) {
        double adjustedCopyNumber = copyNumber.averageTumorCopyNumber();
        builder.adjustedCopyNumber(adjustedCopyNumber);
        double adjustedVAF =
                purityAdjuster.purityAdjustedVAF(copyNumber.chromosome(), Math.max(0.001, adjustedCopyNumber), depth.alleleFrequency());
        builder.adjustedVAF(adjustedVAF);

        double variantPloidy = adjustedCopyNumber * adjustedVAF;
        builder.ploidy(variantPloidy);

        boolean loh = Doubles.lessOrEqual(adjustedCopyNumber, 0) || Doubles.greaterOrEqual(variantPloidy, adjustedCopyNumber - 0.5);
        builder.lossOfHeterozygosity(loh);

        final Clonality clonality = Clonality.fromSample(purityAdjuster, copyNumber, depth);
        builder.clonality(clonality);
    }
}
