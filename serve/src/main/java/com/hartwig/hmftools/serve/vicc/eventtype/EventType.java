package com.hartwig.hmftools.serve.vicc.eventtype;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class EventType {

    public abstract boolean combinedEvent();

    @NotNull
    public abstract String biomarkerType();

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract ViccSource source();

    @NotNull
    public abstract Map<String, List<String>> eventMap();

}
