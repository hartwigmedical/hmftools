package com.hartwig.hmftools.patientreporter.purple;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PurpleAnalysis {

    @NotNull
    public abstract FittedPurity fittedPurity();

    @NotNull
    public abstract FittedPurityScore fittedScorePurity();

    @NotNull
    public abstract List<PurpleCopyNumber> copyNumbers();

    @NotNull
    public List<CopyNumber> ploidyAdjustedCopyNumbers() {
        return copyNumbers().stream().map(x -> ploidyAdjusted(fittedPurity().ploidy(), x)).collect(Collectors.toList());
    }

    @NotNull
    public List<VariantReport> enrich(@NotNull List<VariantReport> variants) {
        return VariantCopyNumberZipper.zip(fittedPurity(), variants, copyNumbers());
    }

    public double purityUncertainty() {
        return fittedScorePurity().maxPurity() - fittedScorePurity().minPurity();
    }

    @NotNull
    @VisibleForTesting
    static CopyNumber ploidyAdjusted(double ploidy, @NotNull PurpleCopyNumber copyNumber) {
        double adjustedCopyNumber = copyNumber.value() <= 1
                ? copyNumber.averageTumorCopyNumber()
                : copyNumber.averageTumorCopyNumber() / ploidy * 2;

        return ImmutablePurpleCopyNumber
                .builder().from(copyNumber).averageTumorCopyNumber(adjustedCopyNumber).build();
    }

}
