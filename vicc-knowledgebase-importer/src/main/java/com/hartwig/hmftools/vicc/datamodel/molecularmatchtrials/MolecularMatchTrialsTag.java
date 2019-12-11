package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsTag {

    @NotNull
    public abstract String facet();

    @NotNull
    public abstract String compositeKey();

    @NotNull
    public abstract String suppress();

    @NotNull
    public abstract String filterType();

    @NotNull
    public abstract String term();

    @NotNull
    public abstract String custom();

    @NotNull
    public abstract String priority();

    @Nullable
    public abstract String alias();

    @Nullable
    public abstract String manualSuppress();

    @Nullable
    public abstract String generatedBy();

    @Nullable
    public abstract String generatedByTerm();

    @Nullable
    public abstract String id();

    @Nullable
    public abstract String manualPriority();
}
