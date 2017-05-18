package com.hartwig.hmftools.common.purple;

import java.io.Serializable;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.freec.FreecCopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class})
public abstract class FittedCopyNumber implements CopyNumber, Serializable {

    public abstract int bafCount();

    public abstract int fittedPloidy();

    public abstract double deviation();

    public abstract double actualBAF();

    public abstract double modelBAF();

    public abstract double bafDeviation();

    public abstract double tumorCNVRatio();

    public abstract double modelCNVRatio();

    public abstract double normalCNVRatio();

    public abstract double cnvDeviation();

    public abstract double normalisedTumorRatio();

    public abstract double ratioOfRatios();

    public abstract double broadRatioOfRatios();

    public abstract double broadBAF();

    public abstract double segmentRatioOfRatios();

    public abstract double segmentBAF();

    public abstract String genotype();

    public abstract FreecCopyNumber.Status status();
}
