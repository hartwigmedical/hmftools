package com.hartwig.hmftools.serve.actionability.fusion;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class ActionableFusionConsolidation {

    private ActionableFusionConsolidation() {
    }

    @NotNull
    public static Set<ActionableFusion> consolidate(@NotNull Set<ActionableFusion> fusions) {
        Map<ActionableFusion, Set<String>> urlsPerEvent = Maps.newHashMap();
        for (ActionableFusion fusion : fusions) {
            ActionableFusion key = stripUrls(fusion);
            Set<String> urls = urlsPerEvent.get(key);
            if (urls == null) {
                urls = Sets.newHashSet();
            }
            urls.addAll(fusion.urls());
            urlsPerEvent.put(key, urls);
        }

        Set<ActionableFusion> consolidatedFusions = Sets.newHashSet();
        for (Map.Entry<ActionableFusion, Set<String>> entry : urlsPerEvent.entrySet()) {
            consolidatedFusions.add(ImmutableActionableFusion.builder().from(entry.getKey()).urls(entry.getValue()).build());
        }
        return consolidatedFusions;
    }

    @NotNull
    private static ActionableFusion stripUrls(@NotNull ActionableFusion fusion) {
        return ImmutableActionableFusion.builder().from(fusion).urls(Sets.newHashSet()).build();
    }
}
