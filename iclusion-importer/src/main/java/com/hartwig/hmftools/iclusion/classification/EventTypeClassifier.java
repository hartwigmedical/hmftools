package com.hartwig.hmftools.iclusion.classification;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.serve.classification.EventMatcher;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.iclusion.data.IclusionMutation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EventTypeClassifier {

    private static final Logger LOGGER = LogManager.getLogger(EventTypeClassifier.class);

    @NotNull
    private static final Map<EventType, EventMatcher> MATCHERS = Maps.newHashMap();

    public EventTypeClassifier() {
    }

    @NotNull
    public static EventType classify(@NotNull IclusionMutation mutation) {
        return EventType.UNKNOWN;
    }
}
