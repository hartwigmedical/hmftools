package com.hartwig.hmftools.serve.actionability.range;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.util.UrlConsolidator;

import org.jetbrains.annotations.NotNull;

public class ActionableRangeUrlConsolidator implements UrlConsolidator<ActionableRange> {

    @NotNull
    @Override
    public ActionableRange stripUrls(@NotNull final ActionableRange instance) {
        return ImmutableActionableRange.builder().from(instance).evidenceUrls(Sets.newHashSet()).build();
    }

    @NotNull
    @Override
    public ActionableRange buildWithUrls(@NotNull final ActionableRange instance, @NotNull final Set<String> urls) {
        return ImmutableActionableRange.builder().from(instance).evidenceUrls(urls).build();
    }
}
