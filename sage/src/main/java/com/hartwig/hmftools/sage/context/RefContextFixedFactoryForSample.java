package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

public class RefContextFixedFactoryForSample {

    private final SageConfig config;
    private final List<VariantHotspot> sortedVariants = Lists.newArrayList();
    private final Map<VariantHotspot, Candidate> readContextMap;

    public RefContextFixedFactoryForSample(@NotNull final SageConfig config, @NotNull final Map<VariantHotspot, Candidate> candidates) {
        this.config = config;
        this.readContextMap = candidates;
        sortedVariants.addAll(candidates.keySet());
        sortedVariants.sort(new VariantHotspotComparator());
    }

    @NotNull
    public List<VariantHotspot> loci() {
        return sortedVariants;
    }

    @NotNull
    public RefContextFixedFactory create(@NotNull final String sample) {
        final RefContextFixedFactory candidates = new RefContextFixedFactory(sample);
        if (sortedVariants.isEmpty()) {
            sortedVariants.addAll(readContextMap.keySet());
            sortedVariants.sort(new VariantHotspotComparator());
        }

        for (VariantHotspot variant : sortedVariants) {
            final Candidate candidate = readContextMap.get(variant);
            candidates.add(variant.chromosome(),
                    variant.position(),
                    variant.ref(),
                    variant.alt(),
                    candidate.maxDepth() < config.maxReadDepth(),
                    candidate.readContext());
        }

        return candidates;
    }
}
