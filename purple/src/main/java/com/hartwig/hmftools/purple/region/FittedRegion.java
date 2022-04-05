package com.hartwig.hmftools.purple.region;

import com.hartwig.hmftools.common.utils.Doubles;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class FittedRegion implements ObservedRegion {

    private static final double MAX_DIPLOID_COPY_NUMBER = 1.2;
    private static final double MIN_DIPLOID_COPY_NUMBER = 0.8;

    public double minorAlleleCopyNumber() {
        return tumorCopyNumber() - majorAlleleCopyNumber();
    }

    public double majorAlleleCopyNumber() {
        return tumorBAF() * tumorCopyNumber();
    }

    public abstract double minorAlleleCopyNumberDeviation();

    public abstract double majorAlleleCopyNumberDeviation();

    public abstract double deviationPenalty();

    public abstract double eventPenalty();

    public abstract double refNormalisedCopyNumber();

    public abstract double tumorCopyNumber();

    public abstract double tumorBAF();

    public abstract double fittedTumorCopyNumber();

    public abstract double fittedBAF();

    public boolean isDiploid() {
        return Doubles.greaterOrEqual(majorAlleleCopyNumber(), MIN_DIPLOID_COPY_NUMBER) && Doubles.lessOrEqual(majorAlleleCopyNumber(),
                MAX_DIPLOID_COPY_NUMBER) && Doubles.greaterOrEqual(minorAlleleCopyNumber(), MIN_DIPLOID_COPY_NUMBER) && Doubles.lessOrEqual(
                minorAlleleCopyNumber(), MAX_DIPLOID_COPY_NUMBER);
    }

}
