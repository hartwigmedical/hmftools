package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EventType {

    @NotNull
    public abstract String gene();

    public abstract boolean combinedEvent();

    public abstract String event();

    @NotNull
    public abstract List<String> interpretEventType();

    @NotNull
    public abstract String biomarkerType();
}
