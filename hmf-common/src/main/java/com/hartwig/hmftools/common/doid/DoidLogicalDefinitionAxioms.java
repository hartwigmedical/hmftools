package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidLogicalDefinitionAxioms
{
    @NotNull
    public abstract String definedClassId();

    @NotNull
    public abstract List<String> genusIds();

    @NotNull
    public abstract List<DoidRestriction> restrictions();
}
