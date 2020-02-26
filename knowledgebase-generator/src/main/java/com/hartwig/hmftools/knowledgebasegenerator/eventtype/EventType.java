package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EventType {

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String eventType();
}
