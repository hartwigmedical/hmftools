package com.hartwig.hmftools.common.refseq;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RefSeq {

    @NotNull
    public abstract String geneId();

    @NotNull
    public abstract String transcriptId();

    @NotNull
    public abstract String displayLabel();

    @NotNull
    public abstract String dbPrimaryAcc();

}
