package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class SequenceOntology {

    @NotNull
    public abstract List<String> hierarchy();

    @NotNull
    public abstract String soid();

    @NotNull
    public abstract String parentSoid();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String parentName();

}
