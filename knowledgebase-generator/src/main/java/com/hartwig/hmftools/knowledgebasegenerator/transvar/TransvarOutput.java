package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TransvarOutput {

    @NotNull
    public abstract String input();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String strand();

    @NotNull
    public abstract String coordinates_gDNA_cDNA_protein();

    @NotNull
    public abstract String region();

    @NotNull
    public abstract String info();

}
