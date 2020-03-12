package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContextFixedFactoryForSample;

import org.jetbrains.annotations.NotNull;

public class ReadContextCollection {

    private final SageConfig config;
    private final Map<VariantHotspot, Candidate> readContextMap = Maps.newHashMap();

    public ReadContextCollection(final SageConfig config) {
        this.config = config;
    }

    public void addCandidates(@NotNull final Collection<AltContext> altContexts) {
        for (final AltContext altContext : altContexts) {
            if (altContext.primaryReadContext().readContext().isComplete()) {
                readContextMap.computeIfAbsent(altContext, x -> new Candidate(altContext)).update(altContext);
            }
        }
    }

    @NotNull
    public RefContextFixedFactoryForSample supplier() {
        return new RefContextFixedFactoryForSample(config, readContextMap);
    }

}
