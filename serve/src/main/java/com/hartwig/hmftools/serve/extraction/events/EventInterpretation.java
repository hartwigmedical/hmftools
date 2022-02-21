package com.hartwig.hmftools.serve.extraction.events;

import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EventInterpretation {

    @NotNull
    public abstract Knowledgebase knowledgebase();

    @NotNull
    public abstract String sourceEvent();

    @NotNull
    public abstract String interpretGene();

    @NotNull
    public abstract String interpretEvent();

    @NotNull
    public abstract EventType interpretEventType();
}