package com.hartwig.hmftools.common.variant.structural.linx;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxCluster
{
    public abstract int clusterId();
    public abstract String resolvedType();
    public abstract boolean synthetic();
    public abstract String subType();
    public abstract int clusterCount();
    public abstract String clusterDesc();
}
