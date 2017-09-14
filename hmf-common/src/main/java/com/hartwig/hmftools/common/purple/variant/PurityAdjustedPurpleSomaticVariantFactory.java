package com.hartwig.hmftools.common.purple.variant;

import static com.hartwig.hmftools.common.purple.variant.ImmutablePurityAdjustedPurpleSomaticVariant.builder;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantBuilder;

import org.jetbrains.annotations.NotNull;

public class PurityAdjustedPurpleSomaticVariantFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;

    public PurityAdjustedPurpleSomaticVariantFactory(@NotNull FittedPurity purity, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        purityAdjuster = new PurityAdjuster(Gender.MALE, purity);
        copyNumberSelector = GenomeRegionSelectorFactory.create(copyNumbers);
    }

    @NotNull
    public List<PurityAdjustedSomaticVariant> create(@NotNull List<PurpleSomaticVariant> variants) {
        final List<PurityAdjustedSomaticVariant> result = Lists.newArrayList();

        for (PurpleSomaticVariant variant : variants) {
            final PurityAdjustedSomaticVariantBuilder builder =
                    builder().from(variant).adjustedCopyNumber(0).adjustedVAF(0).clonality(Clonality.UNKNOWN).lossOfHeterozygosity(false);
            copyNumberSelector.select(variant).ifPresent(x -> builder.purityAdjustment(purityAdjuster, x, variant));
            result.add(builder.build());
        }

        return result;
    }

}
