package com.hartwig.hmftools.serve.actionability.gene;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.util.UrlConsolidator;

import org.jetbrains.annotations.NotNull;

public class ActionableGeneUrlConsolidator implements UrlConsolidator<ActionableGene> {

    @NotNull
    @Override
    public ActionableGene stripUrls(@NotNull final ActionableGene instance) {
        return ImmutableActionableGene.builder().from(instance).evidenceUrls(Sets.newHashSet()).build();
    }

    @NotNull
    @Override
    public ActionableGene buildWithUrls(@NotNull final ActionableGene instance, @NotNull final Set<String> urls) {
        return ImmutableActionableGene.builder().from(instance).evidenceUrls(urls).build();
    }
}
