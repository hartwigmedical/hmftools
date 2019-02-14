package com.hartwig.hmftools.patientreporter.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ReportableGeneDisruption {

    @NotNull
    public abstract String location();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String range();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract Double ploidy();

    @Nullable
    public abstract Integer geneMinCopies();

    @Nullable
    public abstract Integer geneMaxCopies();

    public abstract int firstAffectedExon();
}
