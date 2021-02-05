package com.hartwig.hmftools.ckb.datamodelinterpretation.drug;

import java.util.Date;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Reference {

    public abstract int id();

    @NotNull
    public abstract String pubMedId();

    @NotNull
    public abstract String title();

    @NotNull
    public abstract String url();

    @NotNull
    public abstract String authors();

    @NotNull
    public abstract String journal();

    @NotNull
    public abstract String volume();

    @NotNull
    public abstract String issue();

    @NotNull
    public abstract Date date();

    @NotNull
    public abstract String abstractText();

    @NotNull
    public abstract String year();

}
