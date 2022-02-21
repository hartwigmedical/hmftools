package com.hartwig.hmftools.serve.actionability.util;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.jetbrains.annotations.NotNull;

public final class ActionableEventUrlMerger {

    private ActionableEventUrlMerger() {
    }

    @NotNull
    public static <T extends ActionableEvent> Set<T> merge(@NotNull Set<T> events, @NotNull UrlConsolidator<T> consolidator) {
        Map<T, Set<String>> urlsPerEvent = Maps.newHashMap();
        for (T event : events) {
            T key = consolidator.stripUrls(event);
            Set<String> urls = urlsPerEvent.get(key);
            if (urls == null) {
                urls = Sets.newTreeSet();
            }
            urls.addAll(event.evidenceUrls());
            urlsPerEvent.put(key, urls);
        }

        Set<T> consolidatedEvents = Sets.newHashSet();
        for (Map.Entry<T, Set<String>> entry : urlsPerEvent.entrySet()) {
            consolidatedEvents.add(consolidator.buildWithUrls(entry.getKey(), entry.getValue()));
        }
        return consolidatedEvents;
    }
}
