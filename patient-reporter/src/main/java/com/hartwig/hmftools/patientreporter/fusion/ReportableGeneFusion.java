package com.hartwig.hmftools.patientreporter.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGeneFusion {

    @NotNull
    public abstract String geneStart();

    @NotNull
    public abstract String geneContextStart();

    @NotNull
    public abstract String geneStartTranscript();

    @NotNull
    public abstract String geneEnd();

    @NotNull
    public abstract String geneContextEnd();

    @NotNull
    public abstract String geneEndTranscript();

    public abstract double ploidy();

    @NotNull
    public abstract String source();
}
