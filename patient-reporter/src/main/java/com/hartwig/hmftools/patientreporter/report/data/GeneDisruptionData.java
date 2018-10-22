package com.hartwig.hmftools.patientreporter.report.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class GeneDisruptionData {

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String range();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String copies();

    public abstract int geneMinCopies();

    public abstract int geneMaxCopies();

    public abstract int exonUpstream();
}
