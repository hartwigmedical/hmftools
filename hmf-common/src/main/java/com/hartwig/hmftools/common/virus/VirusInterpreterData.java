package com.hartwig.hmftools.common.virus;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class VirusInterpreterData
{
    @NotNull
    public abstract List<AnnotatedVirus> allViruses();
}
