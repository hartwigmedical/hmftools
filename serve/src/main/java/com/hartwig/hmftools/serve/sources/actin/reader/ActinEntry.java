package com.hartwig.hmftools.serve.sources.actin.reader;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinEventTypeExtractor;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActinEntry {

    @NotNull
    @Value.Derived
    public EventType type() {
        return ActinEventTypeExtractor.extractEventType(this);
    }

    @NotNull
    public abstract String trial();

    @NotNull
    public abstract ActinRule rule();

    @NotNull
    public abstract String gene();

    @Nullable
    public abstract String mutation();
}
