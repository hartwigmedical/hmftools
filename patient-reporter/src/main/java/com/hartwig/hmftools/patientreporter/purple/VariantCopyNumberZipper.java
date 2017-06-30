package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipper;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipperAllPositionsHandler;
import com.hartwig.hmftools.patientreporter.variants.ImmutableVariantReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class VariantCopyNumberZipper implements SimpleGenomeZipperAllPositionsHandler<PurpleCopyNumber, VariantReport> {

    @NotNull
    static List<VariantReport> zip(@NotNull FittedPurity purity, @NotNull final List<VariantReport> variants,
            @NotNull final List<PurpleCopyNumber> copyNumbers) {
        final VariantCopyNumberZipper handler = new VariantCopyNumberZipper(purity);
        SimpleGenomeZipper.zip(copyNumbers, variants, handler);
        return handler.reports();
    }

    private final List<VariantReport> reports = Lists.newArrayList();
    private final PurityAdjuster purityAdjuster;

    private VariantCopyNumberZipper(final FittedPurity purity) {
        this.purityAdjuster = new PurityAdjuster(Gender.MALE, purity);
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

    private VariantReport enrich(@NotNull final PurpleCopyNumber copyNumber, @NotNull final VariantReport variant) {
        double adjustedVAF = Math.min(1, purityAdjuster.purityAdjustedVAF(copyNumber.averageTumorCopyNumber(), variant.alleleFrequency()));
        return ImmutableVariantReport.builder().from(variant).baf(copyNumber.descriptiveBAF()).impliedVAF(adjustedVAF).build();
    }
}
