package com.hartwig.hmftools.common.purple;

import java.io.Serializable;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.freec.FreecStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class})
public abstract class FittedCopyNumber implements CopyNumber, Serializable {

    public abstract int bafCount();

    public abstract int fittedPloidy();

    public abstract double deviation();

    public abstract double observedBAF();

    public abstract double modelBAF();

    public abstract double purityAdjustedBAF();

    public abstract double bafDeviation();

    public abstract double observedNormalRatio();

    public abstract double observedTumorRatio();

    public abstract double modelTumorRatio();

    public abstract double cnvDeviation();

    public abstract double tumorCopyNumber();

    public abstract double broadRatioOfRatios();

    public abstract double broadBAF();

    public abstract double segmentRatioOfRatios();

    public abstract double segmentBAF();

    public abstract FreecStatus status();
}
