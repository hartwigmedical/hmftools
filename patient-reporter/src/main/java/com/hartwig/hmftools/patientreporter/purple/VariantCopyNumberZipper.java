package com.hartwig.hmftools.patientreporter.purple;

import static com.hartwig.hmftools.common.purity.PurityAdjustment.purityAdjustedVAF;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipper;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipperAllPositionsHandler;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class VariantCopyNumberZipper implements SimpleGenomeZipperAllPositionsHandler<PurpleCopyNumber, VariantReport> {

    @NotNull
    static List<VariantReport> zip(@NotNull FittedPurity purity, @NotNull final List<VariantReport> variants, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        VariantCopyNumberZipper handler = new VariantCopyNumberZipper(purity.purity());
        SimpleGenomeZipper.zip(copyNumbers, variants, handler);
        return handler.reports();
    }

    private final List<VariantReport> reports = Lists.newArrayList();
    private final double purity;

    private VariantCopyNumberZipper(final double purity) {
        this.purity = purity;
    }

    private List<VariantReport> reports() {
        return reports;
    }

    @Override
    public void handle(@Nullable final PurpleCopyNumber copyNumber, @NotNull final VariantReport variant) {
        if (copyNumber != null) {
            reports.add(enrich(copyNumber, variant));
        } else {
            reports.add(variant);
        }
    }

    private  VariantReport enrich(@NotNull final PurpleCopyNumber copyNumber,
            @NotNull final VariantReport variant) {
        double adjustedVAF = Math.min(1, purityAdjustedVAF(purity, copyNumber.averageTumorCopyNumber(), variant.alleleFrequency()));
        return ImmutableVariantReport.builder().from(variant).baf(copyNumber.descriptiveBAF()).impliedVAF(adjustedVAF).build();
    }
}
