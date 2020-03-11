package com.hartwig.hmftools.sage.context;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class FixedRefContextCandidatesFactory {

    private final List<VariantHotspot> sortedVariants = Lists.newArrayList();
    private final Map<VariantHotspot, ReadContext> readContextMap;

    public FixedRefContextCandidatesFactory(@NotNull final Map<VariantHotspot, ReadContext> readContextMap) {
        this.readContextMap = readContextMap;
        sortedVariants.addAll(readContextMap.keySet());
        sortedVariants.sort(new VariantHotspotComparator());
    }

    @NotNull
    public List<VariantHotspot> loci() {
        return sortedVariants;
    }

    @NotNull
    public FixedRefContextCandidates create(@NotNull final String sample) {
        final FixedRefContextCandidates candidates = new FixedRefContextCandidates(sample);
        if (sortedVariants.isEmpty()) {
            sortedVariants.addAll(readContextMap.keySet());
            sortedVariants.sort(new VariantHotspotComparator());
        }

        for (VariantHotspot variant : sortedVariants) {
            final ReadContext readContext = readContextMap.get(variant);
            candidates.create(variant.chromosome(), variant.position(), variant.ref(), variant.alt(), readContext);
        }

        return candidates;
    }
}
