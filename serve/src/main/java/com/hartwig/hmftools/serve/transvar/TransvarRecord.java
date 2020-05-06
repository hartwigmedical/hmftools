package com.hartwig.hmftools.serve.transvar;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TransvarRecord {

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String chromosome();

    public abstract long gdnaPosition();

    @NotNull
    public abstract String gdnaRef();

    @NotNull
    public abstract String gdnaAlt();

    // Field is only populated for SNV/MNV
    @Nullable
    public abstract String referenceCodon();

    // Field is only populated for SNV/MNV
    @Nullable
    public abstract List<String> candidateCodons();

    // Field is only populated for Indels
    @Nullable
    public abstract Integer indelLength();

}
