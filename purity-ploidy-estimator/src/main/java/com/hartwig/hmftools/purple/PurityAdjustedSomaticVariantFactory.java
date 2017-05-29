package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedSomaticVariants;

import java.util.List;
import java.util.function.BiConsumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.ConsolidatedRegion;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.zipper.GenomeZipper;
import com.hartwig.hmftools.common.zipper.GenomeZipperRegionHandler;

class PurityAdjustedSomaticVariantFactory implements GenomeZipperRegionHandler<ConsolidatedRegion>, BiConsumer<ConsolidatedRegion, SomaticVariant> {

    static List<SomaticVariant> purpleAdjusted(double purity, List<ConsolidatedRegion> regions, List<SomaticVariant> variants) {

        PurityAdjustedSomaticVariantFactory zipHandler = new PurityAdjustedSomaticVariantFactory(purity);
        GenomeZipper<ConsolidatedRegion> zipper = new GenomeZipper<>(regions, zipHandler);
        zipper.addPositions(variants, zipHandler);
        zipper.run();
        return zipHandler.results();
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
    public void enter(final ConsolidatedRegion region) {
        // Ignore
    }

    @Override
    public void exit(final ConsolidatedRegion region) {
        // Ignore
    }

    @Override
    public void accept(final ConsolidatedRegion consolidatedRegion, final SomaticVariant variant) {
        if (variant.type() == VariantType.SNP) {
            results.add(purityAdjusted(consolidatedRegion.averageTumorCopyNumber(), variant));
        }
    }

    private SomaticVariant purityAdjusted(double ploidy, SomaticVariant variant) {
        double purityAdjustedFrequency = purityAdjustedSomaticVariants(purity, ploidy, variant.alleleFrequency());
        return SomaticVariant.Builder.fromVariant(variant)
                .alleleReadCount((int) Math.round(purityAdjustedFrequency * 10000)).totalReadCount(10000).build();

    }
}
