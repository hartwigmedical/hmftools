package com.hartwig.hmftools.ckb.datamodelinterpretation.reference;

import java.util.Date;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Reference {

    public abstract int id();

    @Nullable
    public abstract String pubMedId();

    @Nullable
    public abstract String title();

    @Nullable
    public abstract String url();

    @Nullable
    public abstract String authors();

    @Nullable
    public abstract String journal();

    @Nullable
    public abstract String volume();

    @Nullable
    public abstract String issue();

    @Nullable
    public abstract String date();

    @Nullable
    public abstract String abstractText();

    @Nullable
    public abstract String year();
}