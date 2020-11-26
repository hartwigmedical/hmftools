package com.hartwig.hmftools.common.serve.classification;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.matchers.EventMatcher;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EventClassifier {

    private static final Logger LOGGER = LogManager.getLogger(EventClassifier.class);

    @NotNull
    private final Map<MutationType, EventMatcher> matchers;

    public EventClassifier(@NotNull final Map<MutationType, EventMatcher> matchers) {
        this.matchers = matchers;
    }

    @NotNull
    public MutationType determineType(@NotNull String gene, @NotNull String event) {
        Map<MutationType, Boolean> evaluations = Maps.newHashMap();

        for (Map.Entry<MutationType, EventMatcher> entry : matchers.entrySet()){
            evaluations.put(entry.getKey(), entry.getValue().matches(gene, event));
        }

        Set<MutationType> positiveTypes = Sets.newHashSet();
        for (Map.Entry<MutationType, Boolean> evaluation : evaluations.entrySet()) {
            if (evaluation.getValue()) {
                positiveTypes.add(evaluation.getKey());
            }
        }

        if (positiveTypes.size() > 1) {
            LOGGER.warn("More than one type evaluated to true for '{}' on '{}': {}", event, gene, positiveTypes);
        } else if (positiveTypes.size() == 1) {
            return positiveTypes.iterator().next();
        }

        return MutationType.UNKNOWN;
    }
}
