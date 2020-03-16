package com.hartwig.hmftools.knowledgebasegenerator.eventtype;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;

import org.apache.logging.log4j.util.Strings;
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
    public abstract Source source();

    @NotNull
    public abstract Map<String, List<String>> eventMap();

    public abstract String event();

    @NotNull
    public abstract List<String> interpretEventType();

}
