package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidEntry
{
    @NotNull
    public abstract String id();

    @NotNull
    public abstract List<DoidNode> nodes();

    @NotNull
    public abstract List<DoidEdge> edges();

    @NotNull
    public abstract DoidGraphMetaData meta();

    @Nullable
    public abstract List<DoidLogicalDefinitionAxioms> logicalDefinitionAxioms();

}
