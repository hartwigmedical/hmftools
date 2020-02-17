package com.hartwig.hmftools.knowledgebasegenerator.transvar;

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

    public abstract long dDNAPosition();

    @NotNull
    public abstract String gDNARef();

    @NotNull
    public abstract String gDNAAlt();

    @NotNull
    public abstract String referenceCodon();

    @NotNull
    public abstract List<String> candidateCodons();


}
