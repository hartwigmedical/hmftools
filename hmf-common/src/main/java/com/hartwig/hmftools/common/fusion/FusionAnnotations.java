package com.hartwig.hmftools.common.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
public abstract class FusionAnnotations
{
    public abstract int clusterId();
    public abstract int clusterCount();
    public abstract String resolvedType();

    // chained fusion info
    @Nullable
    public abstract FusionChainInfo chainInfo();

    // info about whether the fusion is valid to the end of the transcript via the chain
    @Nullable
    public abstract FusionTermination disruptionUp();

    @Nullable
    public abstract FusionTermination disruptionDown();

}
