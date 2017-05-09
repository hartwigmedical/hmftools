package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import org.immutables.value.Value;

import java.io.Serializable;

@Value.Style(allParameters = true)
@Value.Immutable
public abstract class FittedCopyNumber implements CopyNumber, Serializable {

    public abstract int fittedPloidy();

    public abstract double deviation();

    public abstract double actualBAF();

    public abstract double modelBAF();

    public abstract double bafDeviation();

    public abstract double actualCNVRatio();

    public abstract double modelCNVRatio();

    public abstract double cnvDeviation();
}
