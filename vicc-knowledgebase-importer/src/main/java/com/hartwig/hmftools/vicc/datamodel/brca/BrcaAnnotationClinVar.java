package com.hartwig.hmftools.vicc.datamodel.brca;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class BrcaAnnotationClinVar {

    @NotNull
    public abstract String variantInClinVar();

    @NotNull
    public abstract String bxId();

    @NotNull
    public abstract String clinicalSignificance();

    @NotNull
    public abstract String submitter();

    @NotNull
    public abstract String method();

    @NotNull
    public abstract String alleleOrigin();

    @NotNull
    public abstract String scv();

    @NotNull
    public abstract String dateLastUpdated();
}
