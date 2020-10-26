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

    @Nullable
    public abstract String idNodes();

    @Nullable
    public abstract List<String> metaNodes();

    @Nullable
    public abstract List<DoidEquivalentNodesSets> equivalentNodesSets();

    @Nullable
    public abstract List<DoidLogicalDefinitionAxioms> logicalDefinitionAxioms();

    @Nullable
    public abstract List<DoidDomainRangeAxioms> domainRangeAxioms();

    @Nullable
    public abstract List<DoidPropertyChainAxioms> propertyChainAxioms();
}
