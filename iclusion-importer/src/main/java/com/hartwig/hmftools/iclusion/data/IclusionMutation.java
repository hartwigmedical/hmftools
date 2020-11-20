package com.hartwig.hmftools.iclusion.data;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.iclusion.classification.EventTypeClassifier;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class IclusionMutation {

    @NotNull
    @Value.Derived
    public EventType type() {
        return EventTypeClassifier.classify(this);
    }

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String name();

    public abstract boolean negation();
}
