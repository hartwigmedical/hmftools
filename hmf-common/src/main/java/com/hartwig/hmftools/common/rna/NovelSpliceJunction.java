package com.hartwig.hmftools.common.rna;

import org.immutables.value.Value;

@Value.Immutable
public abstract class NovelSpliceJunction
{
    public abstract String geneName();
    public abstract String chromosome();
    public abstract int junctionStart();
    public abstract int junctionEnd();
    public abstract AltSpliceJunctionType type();
    public abstract String transcriptStart();
    public abstract String transcriptEnd();
    public abstract int exonStart();
    public abstract int exonEnd();
    public abstract int fragmentCount();
    public abstract int depthStart();
    public abstract int depthEnd();
    public abstract AltSpliceJunctionContext regionStart();
    public abstract AltSpliceJunctionContext regionEnd();
    public abstract String basesStart();
    public abstract String basesEnd();
    public abstract int cohortFrequency();
}
