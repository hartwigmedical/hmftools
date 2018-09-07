package com.hartwig.hmftools.common.purple.region;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class})
public abstract class FittedRegion implements ObservedRegion {

    public abstract int modelPloidy();

    public abstract double modelBAF();

    public abstract double modelTumorRatio();

    public abstract double deviation();

    public abstract double ploidyPenalty();

    public abstract double bafDeviation();

    public abstract double cnvDeviation();

    public abstract double refNormalisedCopyNumber();

    public abstract double tumorCopyNumber();

    public abstract double tumorBAF();

    public abstract double fittedTumorCopyNumber();

    public abstract double fittedBAF();

    public double majorAllelePloidy() {
        return tumorBAF() * tumorCopyNumber();
    }

    public double minorAllelePloidy() {
        return tumorCopyNumber() - majorAllelePloidy();
    }
}
