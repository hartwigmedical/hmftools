package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedFrequency;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipper;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipperInRegionPositionsHandler;

import org.jetbrains.annotations.NotNull;

class PurityAdjustedSomaticVariantFactory
        implements SimpleGenomeZipperInRegionPositionsHandler<PurpleCopyNumber, SomaticVariant> {

    static List<SomaticVariant> purpleAdjusted(double purity, List<PurpleCopyNumber> regions,
            List<SomaticVariant> variants) {
        PurityAdjustedSomaticVariantFactory handler = new PurityAdjustedSomaticVariantFactory(purity);
        SimpleGenomeZipper.zip(regions, variants, handler);
        return handler.results();
    }

    private final double purity;
    private final List<SomaticVariant> results = Lists.newArrayList();

    private PurityAdjustedSomaticVariantFactory(final double purity) {
        this.purity = purity;
    }

    private List<SomaticVariant> results() {
        return results;
    }

    @Override
    public void handle(@NotNull final PurpleCopyNumber consolidatedRegion, @NotNull final SomaticVariant variant) {
        if (variant.type() == VariantType.SNP) {
            results.add(purityAdjusted(consolidatedRegion.averageTumorCopyNumber(), variant));
        }
    }

    private SomaticVariant purityAdjusted(double ploidy, SomaticVariant variant) {
        double purityAdjustedFrequency = purityAdjustedSomaticVariants(purity, ploidy, variant.alleleFrequency());
        return SomaticVariant.Builder.fromVariant(variant)
                .alleleReadCount((int) Math.round(purityAdjustedFrequency * 10000))
                .totalReadCount(10000)
                .build();
    }

    private static double purityAdjustedSomaticVariants(final double purity, final double ploidy,
            final double observedVAF) {
        return purityAdjustedFrequency(purity, ploidy, observedVAF, 0);
    }
}
