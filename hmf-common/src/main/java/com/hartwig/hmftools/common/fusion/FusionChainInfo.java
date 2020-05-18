package com.hartwig.hmftools.common.fusion;

import org.immutables.value.Value;

@Value.Immutable
public abstract class FusionChainInfo
{
    public abstract int chainId();
    public abstract int chainLinks();
    public abstract int chainLength();
    public abstract boolean validTraversal();
    public abstract boolean traversalAssembled();
}
