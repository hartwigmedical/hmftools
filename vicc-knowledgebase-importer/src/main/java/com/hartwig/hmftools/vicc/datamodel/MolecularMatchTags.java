package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTags {

    @NotNull
    public abstract String priority();

    @Nullable
    public abstract String compositeKey();

    @Nullable
    public abstract String suppress();

    @Nullable
    public abstract String filterType();

    @NotNull
    public abstract String term();

    @Nullable
    public abstract String primary();

    @NotNull
    public abstract String facet();

    @Nullable
    public abstract String valid();

    @Nullable
    public abstract String custom();

    @Nullable
    public abstract String isNew();

    @Nullable
    public abstract String generatedBy();

    @Nullable
    public abstract String manualSuppress();

    @Nullable
    public abstract String generatedByTerm();

    @Nullable
    public abstract String transcript();
}
