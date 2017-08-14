package com.hartwig.hmftools.common.purple.copynumber;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;

final class CopyNumberBuilder {

    private static final double MIN_COPY_NUMBER_TOLERANCE = 0.3;
    @VisibleForTesting
    static final double LC_MAX_COPY_NUMBER_TOLERANCE = 1.3;
    @VisibleForTesting
    static final double HC_MAX_COPY_NUMBER_TOLERANCE = 0.7;

    @NotNull
    private final CombinedFittedRegion combinedRegion;
    private final CopyNumberDeviation deviation;

    CopyNumberBuilder(boolean highConfidence, @NotNull final PurityAdjuster purityAdjuster, final FittedRegion combinedRegion) {
        this.combinedRegion = new CombinedFittedRegion(highConfidence, combinedRegion);
        this.deviation = new CopyNumberDeviation(purityAdjuster);
    }

    public String chromosome() {
        return combinedRegion.region().chromosome();
    }

    public final double averageObservedBAF() {
        return combinedRegion.region().observedBAF();
    }

    public final double averageTumorCopyNumber() {
        return combinedRegion.region().tumorCopyNumber();
    }


    public void extendRegion(@NotNull final FittedRegion value) {
        assert (chromosome().equals(value.chromosome())) : "Regions cannot be extended between chromosomes";
        combinedRegion.combine(value);
    }

    @NotNull
    public FittedRegion build() {
        return combinedRegion.region();
    }

    boolean withinCopyNumberTolerance(@NotNull final FittedRegion copyNumber) {
        return deviation.withinCopyNumberTolerance(build(), copyNumber);
    }

}
