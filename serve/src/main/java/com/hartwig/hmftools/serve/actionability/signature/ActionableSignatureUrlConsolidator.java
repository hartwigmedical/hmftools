package com.hartwig.hmftools.serve.actionability.signature;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.util.UrlConsolidator;

import org.jetbrains.annotations.NotNull;

public class ActionableSignatureUrlConsolidator implements UrlConsolidator<ActionableSignature> {

    @NotNull
    @Override
    public ActionableSignature stripUrls(@NotNull final ActionableSignature instance) {
        return ImmutableActionableSignature.builder().from(instance).urls(Sets.newHashSet()).build();
    }

    @NotNull
    @Override
    public ActionableSignature buildWithUrls(@NotNull final ActionableSignature instance, @NotNull final Set<String> urls) {
        return ImmutableActionableSignature.builder().from(instance).urls(urls).build();
    }
}
