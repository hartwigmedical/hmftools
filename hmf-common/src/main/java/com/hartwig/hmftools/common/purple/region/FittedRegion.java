package com.hartwig.hmftools.common.purple.region;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class})
public abstract class FittedRegion implements ObservedRegion {

    public abstract int fittedPloidy();

    public abstract double deviation();

    public abstract double modelBAF();

    public abstract double bafDeviation();

    public abstract double modelTumorRatio();

    public abstract double cnvDeviation();

    public abstract double tumorCopyNumber();

    public abstract double refNormalisedCopyNumber();

    public abstract double broadTumorCopyNumber();

    public abstract double broadBAF();

    public abstract double segmentTumorCopyNumber();

    public abstract double segmentBAF();
}
