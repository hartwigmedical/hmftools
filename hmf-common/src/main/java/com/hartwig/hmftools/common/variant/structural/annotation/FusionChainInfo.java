package com.hartwig.hmftools.common.variant.structural.annotation;

import org.immutables.value.Value;

@Value.Immutable
public abstract class FusionChainInfo
{
    public abstract int chainId();
    public abstract int chainLinks();
    public abstract long chainLength();
    public abstract boolean validTraversal();
    public abstract boolean traversalAssembled();
}
