package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidEntry {

    @NotNull
    public abstract List<DoidNode> doidNodes();

    @NotNull
    public abstract DoidEdges edges();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract DoidGraphMetaData meta();

    @NotNull
    public abstract List<DoidEquivalentNodesSets> equivalentNodesSets();

    @NotNull
    public abstract List<DoidLogicalDefinitionAxioms> logicalDefinitionAxioms();

    @NotNull
    public abstract List<DoidDomainRangeAxioms> domainRangeAxioms();

    @NotNull
    public abstract List<DoidPropertyChainAxioms> propertyChainAxioms();

}
