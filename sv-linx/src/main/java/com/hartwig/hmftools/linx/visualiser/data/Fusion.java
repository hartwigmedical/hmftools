package com.hartwig.hmftools.linx.visualiser.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Fusion
{
    public abstract String sampleId();

    public abstract int clusterId();

    public abstract String geneUp();

    public abstract String chromosomeUp();

    public abstract long positionUp();

    public abstract String regionTypeUp();

    public abstract int strandUp();

    public abstract int exonUp();

    public abstract String geneDown();

    public abstract String chromosomeDown();

    public abstract long positionDown();

    public abstract String regionTypeDown();

    public abstract int strandDown();

    public abstract int exonDown();

}
