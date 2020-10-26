package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidNodes {

    @Nullable
    public abstract List<DoidEdge> edges();

    @NotNull
    public abstract String idNodes();

    @NotNull
    public abstract List<String> metaNodes();

    @NotNull
    public abstract List<String> equivalentNodesSets();

    @NotNull
    public abstract List<String> logicalDefinitionAxioms();

    @NotNull
    public abstract List<String> domainRangeAxioms();

    @NotNull
    public abstract List<String> propertyChainAxioms();
}
