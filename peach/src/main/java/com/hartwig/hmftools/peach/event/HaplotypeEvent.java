package com.hartwig.hmftools.peach.event;

import org.jetbrains.annotations.NotNull;

public interface HaplotypeEvent
{
    String EVENT_ID_DELIMITER = "_";

    @NotNull
    String id();

    boolean isRelevantFor(@NotNull HaplotypeEvent event);
}
