package com.hartwig.hmftools.common.purple.region;

import com.hartwig.hmftools.common.numeric.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class FittedRegion implements ObservedRegion {

    private static final double MAX_DIPLOID_PLOIDY = 1.2;
    private static final double MIN_DIPLOID_PLOIDY = 0.8;

    public abstract double minorAllelePloidy();

    public abstract double majorAllelePloidy();

    public abstract double minorAllelePloidyDeviation();

    public abstract double majorAllelePloidyDeviation();

    public abstract double deviation();

    public abstract double ploidyPenalty();

    public abstract double refNormalisedCopyNumber();

    public abstract double tumorCopyNumber();

    public abstract double tumorBAF();

    public abstract double fittedTumorCopyNumber();

    public abstract double fittedBAF();

    public boolean isDiploid() {
        return Doubles.greaterOrEqual(majorAllelePloidy(), MIN_DIPLOID_PLOIDY) &&
                Doubles.lessOrEqual(majorAllelePloidy(), MAX_DIPLOID_PLOIDY) &&
                Doubles.greaterOrEqual(minorAllelePloidy(), MIN_DIPLOID_PLOIDY) &&
                Doubles.lessOrEqual(minorAllelePloidy(), MAX_DIPLOID_PLOIDY);

    }

}
