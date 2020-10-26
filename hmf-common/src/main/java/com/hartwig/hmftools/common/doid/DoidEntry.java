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
    public abstract String doid();

    @NotNull
    public abstract String url();

    @Nullable
    public abstract String doidTerm();

    @Nullable
    public abstract String type();

    @Nullable
    public abstract DoidMetadata doidMetadata();

    @NotNull
    public abstract List<String> edges();

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
