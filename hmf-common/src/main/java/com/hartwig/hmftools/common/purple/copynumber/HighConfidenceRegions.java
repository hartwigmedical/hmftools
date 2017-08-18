package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class HighConfidenceRegions {

    private static final int MIN_BAF_COUNT = 50;
    private static final double MAX_BAF_DEVIATION = 0.03;
    private static final double MAX_COPY_NUMBER_DEVIATION = 0.30;

    @NotNull
    private final List<FittedRegion> result = Lists.newArrayList();
    @Nullable
    private CombinedFittedRegion builder;
    @Nullable
    private FittedRegion last;
    private final double maxDeviation;

    HighConfidenceRegions(@NotNull final PurityAdjuster purityAdjuster) {
        maxDeviation = purityAdjuster.purityAdjustedMaxCopyNumberDeviation(MAX_COPY_NUMBER_DEVIATION);
    }

    @NotNull
    List<FittedRegion> highConfidence(@NotNull final List<FittedRegion> fittedRegions) {
        fittedRegions.stream()
                .filter(copyNumber -> copyNumber.bafCount() > MIN_BAF_COUNT && copyNumber.status() == FreecStatus.SOMATIC
                        && !Doubles.isZero(copyNumber.observedTumorRatio()))
                .forEach(this::process);

        endRegion();
        return result;
    }

    private void process(@NotNull final FittedRegion current) {
        if (builder == null || isNewChromosome(current, last) || isLargeDeviation(current)) {
            endRegion();
            builder = new CombinedFittedRegion(true, current);
        } else {
            assert builder != null;
            builder.combine(current);
        }

        last = current;
    }

    private void endRegion() {
        if (builder != null) {
            result.add(builder.region());
            builder = null;
        }
    }

    private static boolean isNewChromosome(@NotNull final FittedRegion current, @Nullable final FittedRegion previous) {
        return previous != null && !current.chromosome().equals(previous.chromosome());
    }

    private boolean isLargeDeviation(@NotNull final FittedRegion current) {
        assert builder != null;

        double copyNumberDeviation = Math.abs(current.tumorCopyNumber() - builder.region().tumorCopyNumber());
        if (Doubles.greaterThan(copyNumberDeviation, maxDeviation)) {
            return true;
        }

        double bafDeviation = Math.abs(current.observedBAF() - builder.region().observedBAF());
        return Doubles.greaterThan(bafDeviation, MAX_BAF_DEVIATION);
    }
}
