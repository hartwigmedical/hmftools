package com.hartwig.hmftools.linx.visualiser.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Fusion
{
    @NotNull
    public abstract String sampleId();

    public abstract int clusterId();

    @NotNull
    public abstract String geneUp();

    @NotNull
    public abstract String chromosomeUp();

    public abstract long positionUp();

    @NotNull
    public abstract String regionTypeUp();

    public abstract int strandUp();

    public abstract int exonUp();

    @NotNull
    public abstract String geneDown();

    @NotNull
    public abstract String chromosomeDown();

    public abstract long positionDown();

    @NotNull
    public abstract String regionTypeDown();

    public abstract int strandDown();

    public abstract int exonDown();

    @NotNull
    public String name() {
        return geneUp() + "_" + geneDown();
    }

}
