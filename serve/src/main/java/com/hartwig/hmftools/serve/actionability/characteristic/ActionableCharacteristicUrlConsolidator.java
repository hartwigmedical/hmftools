package com.hartwig.hmftools.serve.actionability.characteristic;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.actionability.util.UrlConsolidator;

import org.jetbrains.annotations.NotNull;

public class ActionableCharacteristicUrlConsolidator implements UrlConsolidator<ActionableCharacteristic> {

    @NotNull
    @Override
    public ActionableCharacteristic stripUrls(@NotNull final ActionableCharacteristic instance) {
        return ImmutableActionableCharacteristic.builder().from(instance).evidenceUrls(Sets.newHashSet()).build();
    }

    @NotNull
    @Override
    public ActionableCharacteristic buildWithUrls(@NotNull final ActionableCharacteristic instance, @NotNull final Set<String> urls) {
        return ImmutableActionableCharacteristic.builder().from(instance).evidenceUrls(urls).build();
    }
}
