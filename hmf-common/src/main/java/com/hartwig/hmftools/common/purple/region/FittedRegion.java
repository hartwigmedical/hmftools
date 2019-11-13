package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.math.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class FittedRegion implements ObservedRegion {

    private static final double MAX_DIPLOID_PLOIDY = 1.2;
    private static final double MIN_DIPLOID_PLOIDY = 0.8;

    public double minorAllelePloidy() {
        return tumorCopyNumber() - majorAllelePloidy();
    }

    public double majorAllelePloidy() {
        return tumorBAF() * tumorCopyNumber();
    }

    public abstract double minorAllelePloidyDeviation();

    public abstract double majorAllelePloidyDeviation();

    public abstract double deviationPenalty();

    public abstract double eventPenalty();

    public abstract double refNormalisedCopyNumber();

    public abstract double tumorCopyNumber();

    public abstract double tumorBAF();

    public abstract double fittedTumorCopyNumber();

    public abstract double fittedBAF();

    public boolean isDiploid() {
        return Doubles.greaterOrEqual(majorAllelePloidy(), MIN_DIPLOID_PLOIDY) && Doubles.lessOrEqual(majorAllelePloidy(),
                MAX_DIPLOID_PLOIDY) && Doubles.greaterOrEqual(minorAllelePloidy(), MIN_DIPLOID_PLOIDY) && Doubles.lessOrEqual(
                minorAllelePloidy(),
                MAX_DIPLOID_PLOIDY);
    }

}
