package com.hartwig.hmftools.common.purple.region;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.freec.FreecStatus;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.FittedCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class HighConfidenceRegions {

    private static final int MIN_BAF_COUNT = 50;
    private static final double MAX_BAF_DEVIATION = 0.03;
    private static final double MAX_COPY_NUMBER_DEVIATION = 0.30;

    @NotNull
    private final List<ConsolidatedRegion> result = Lists.newArrayList();

    @Nullable
    private ConsolidatedRegionBuilder builder;
    @Nullable
    private FittedCopyNumber last;

    List<ConsolidatedRegion> highConfidence(List<FittedCopyNumber> copyNumbers) {
        for (FittedCopyNumber copyNumber : copyNumbers) {
            if (copyNumber.bafCount() > MIN_BAF_COUNT && copyNumber.status() == FreecStatus.SOMATIC && !Doubles.isZero(
                    copyNumber.observedTumorRatio())) {
                process(copyNumber);
            }
        }

        endRegion();
        return result;
    }

    private void process(@NotNull FittedCopyNumber current) {
        if (builder == null || isNewChromosome(current, last) || isLargeDeviation(current)) {
            endRegion();
            builder = new ConsolidatedRegionBuilder(current);
        } else {
            assert builder != null;
            builder.extendRegion(current);
        }

        last = current;
    }

    private void endRegion() {
        if (builder != null) {
            result.add(builder.build());
            builder = null;
        }
    }

    private boolean isNewChromosome(@NotNull FittedCopyNumber current, @Nullable FittedCopyNumber previous) {
        return previous != null && !current.chromosome().equals(previous.chromosome());
    }

    private boolean isLargeDeviation(@NotNull FittedCopyNumber current) {
        assert builder != null;

        double copyNumberDeviation = Math.abs(current.tumorCopyNumber() - builder.averageTumorCopyNumber());
        if (Doubles.greaterThan(copyNumberDeviation, MAX_COPY_NUMBER_DEVIATION)) {
            return true;
        }

        double bafDeviation = Math.abs(current.observedBAF() - builder.averageObservedBAF());
        return Doubles.greaterThan(bafDeviation, MAX_BAF_DEVIATION);
    }
}
