package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipper;
import com.hartwig.hmftools.common.zipper.SimpleGenomeZipperAllPositionsHandler;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class VariantCopyNumberZipper implements SimpleGenomeZipperAllPositionsHandler<PurpleCopyNumber, VariantReport> {

    @NotNull
    public static List<VariantReport> zip(@NotNull final List<VariantReport> variants, @NotNull final List<PurpleCopyNumber> copyNumbers) {
        VariantCopyNumberZipper handler = new VariantCopyNumberZipper();
        SimpleGenomeZipper.zip(copyNumbers, variants, handler);
        return handler.reports();
    }

    private final List<VariantReport> reports = Lists.newArrayList();

    private VariantCopyNumberZipper() {
        // Empty
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

    private static VariantReport enrich(@NotNull final PurpleCopyNumber copyNumber,
            @NotNull final VariantReport variant) {
        return ImmutableVariantReport.builder().from(variant).baf(copyNumber.descriptiveBAF()).build();
    }
}
