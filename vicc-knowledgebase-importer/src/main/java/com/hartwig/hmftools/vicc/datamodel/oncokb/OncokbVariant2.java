package com.hartwig.hmftools.vicc.datamodel.oncokb;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OncokbVariant2 {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String alteration();

    @NotNull
    public abstract OncoKbConsequence consequence();

    @NotNull
    public abstract OncokbGene2 gene();

    @NotNull
    public abstract String proteinStart();

    @NotNull
    public abstract String proteinEnd();

    @Nullable
    public abstract String refResidues();

    @Nullable
    public abstract String variantResidues();
}
