package com.hartwig.hmftools.common.serve.classification;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ComprehensiveEventClassifier {

    private static final Logger LOGGER = LogManager.getLogger(ComprehensiveEventClassifier.class);

    @NotNull
    private final Map<EventType, EventMatcher> matchers;

    public ComprehensiveEventClassifier(@NotNull final Map<EventType, EventMatcher> matchers) {
        this.matchers = matchers;
    }

    @NotNull
    public EventType determineType(@NotNull String gene, @NotNull String event) {
        Map<EventType, Boolean> evaluations = Maps.newHashMap();

        for (Map.Entry<EventType, EventMatcher> entry : matchers.entrySet()){
            evaluations.put(entry.getKey(), entry.getValue().matches(gene, event));
        }

        Set<EventType> positiveTypes = Sets.newHashSet();
        for (Map.Entry<EventType, Boolean> evaluation : evaluations.entrySet()) {
            if (evaluation.getValue()) {
                positiveTypes.add(evaluation.getKey());
            }
        }

        if (positiveTypes.size() > 1) {
            LOGGER.warn("More than one type evaluated to true for '{}' on '{}': {}", event, gene, positiveTypes);
        } else if (positiveTypes.size() == 1) {
            return positiveTypes.iterator().next();
        }

        return EventType.UNKNOWN;
    }
}
