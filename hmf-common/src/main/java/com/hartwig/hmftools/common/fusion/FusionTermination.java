package com.hartwig.hmftools.common.fusion;

import org.immutables.value.Value;

@Value.Immutable
public abstract class FusionTermination
{
    public abstract int totalBreakends();
    public abstract int facingBreakends();
    public abstract int disruptedExons();
    public abstract boolean transcriptTerminated();
    public abstract long minDistance();
    public abstract boolean allLinksAssembled();

}
