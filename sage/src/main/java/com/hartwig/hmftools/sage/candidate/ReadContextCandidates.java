package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContextFixedFactorySupplier;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

public class ReadContextCandidates {

    private final Map<VariantHotspot, ReadContextCounter> readContextMap = Maps.newHashMap();

    public void addCandidates(@NotNull final Collection<AltContext> altContexts) {
        for (final AltContext altContext : altContexts) {
            readContextMap.merge(altContext, altContext.primaryReadContext(), this::mostSupport);
        }
    }

    @NotNull
    public RefContextFixedFactorySupplier candidateFactory() {
        return new RefContextFixedFactorySupplier(candidateMap());
    }

    @NotNull
    private ReadContextCounter mostSupport(@NotNull final ReadContextCounter oldValue, @NotNull final ReadContextCounter newValue) {
        return oldValue.altSupport() > newValue.altSupport() ? oldValue : newValue;
    }

    @NotNull
    private Map<VariantHotspot, ReadContext> candidateMap() {
        return readContextMap.entrySet().stream().collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().readContext()));
    }
}
