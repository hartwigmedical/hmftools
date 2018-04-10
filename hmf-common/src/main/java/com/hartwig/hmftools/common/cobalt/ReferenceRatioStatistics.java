package com.hartwig.hmftools.common.cobalt;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface ReferenceRatioStatistics {

    double TWO_X_CUTOFF = 0.65;

    double Y_MEDIAN_CUTOFF = 0.05;
    int Y_COUNT_CUTOFF = 1000;


    int xCount();

    double xMedian();

    int yCount();

    double yMedian();

    default boolean containsYChromosome() {
        return Doubles.greaterThan(yMedian(), Y_MEDIAN_CUTOFF) && yCount() > Y_COUNT_CUTOFF;
    }

    default boolean containsTwoXChromosomes() {
        return Doubles.greaterThan(xMedian(), TWO_X_CUTOFF);
    }

}