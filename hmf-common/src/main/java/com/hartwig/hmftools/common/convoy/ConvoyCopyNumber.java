package com.hartwig.hmftools.common.convoy;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import org.immutables.value.Value;

@Value.Style(allParameters = true)
@Value.Immutable
public abstract class ConvoyCopyNumber implements CopyNumber {

    public abstract double mBAF();

    public abstract int mBAFCount();

    public abstract double tumorRatio();

    public abstract double ratioOfRatio();
}
