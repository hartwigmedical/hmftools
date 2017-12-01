package com.hartwig.hmftools.svannotation.annotations;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class GeneDisruption {
    public abstract String geneName();

    public abstract String location();

    public abstract String geneContext();

    public abstract String transcript();

    public abstract String partner();

    public abstract String type();

    public abstract String orientation();

    public abstract String vaf();
}
