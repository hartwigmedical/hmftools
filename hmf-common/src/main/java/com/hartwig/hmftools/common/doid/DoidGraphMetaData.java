package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidGraphMetaData
{
    @Nullable
    public abstract List<String> subsets();

    @Nullable
    public abstract List<DoidXref> xrefs();

    @Nullable
    public abstract List<DoidBasicPropertyValue> basicPropertyValues();

    @Nullable
    public abstract String version();
}
