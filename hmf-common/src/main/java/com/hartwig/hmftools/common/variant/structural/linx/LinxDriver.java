package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxDriver
{
    public abstract int clusterId();
    public abstract String gene();
    public abstract String eventType();

}
