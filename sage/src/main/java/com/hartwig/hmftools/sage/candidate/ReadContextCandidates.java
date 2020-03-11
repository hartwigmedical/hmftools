package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.FixedRefContextCandidatesFactory;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;

import org.jetbrains.annotations.NotNull;

public class ReadContextCandidates {

    private final Map<VariantHotspot, ReadContextCounter> readContextMap = Maps.newHashMap();

    public void addCandidates(@NotNull final Collection<AltContext> altContexts) {
        for (final AltContext altContext : altContexts) {
            final VariantHotspot variant = ImmutableVariantHotspotImpl.builder().from(altContext).build();
            readContextMap.merge(variant, altContext.primaryReadContext(), this::mostSupport);
        }
    }


    @NotNull
    public FixedRefContextCandidatesFactory candidateFactory() {
        return new FixedRefContextCandidatesFactory(candidateMap());
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
