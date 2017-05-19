package com.hartwig.hmftools.common.purple;

import java.io.Serializable;

import com.hartwig.hmftools.common.copynumber.CopyNumber;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true, passAnnotations = { NotNull.class, Nullable.class})
public abstract class EnrichedCopyNumber implements CopyNumber, Serializable {

    public abstract double mBAF();

    public abstract int mBAFCount();

    public abstract double tumorRatio();

    public abstract double normalRatio();

    public abstract double ratioOfRatio();

}
