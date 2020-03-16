package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.context.AltContext;

import org.jetbrains.annotations.NotNull;

public class Candidates {

    private final Map<VariantHotspot, Candidate> candidateMap = Maps.newHashMap();
    private List<Candidate> candidateList;

    public void add(@NotNull final Collection<AltContext> altContexts) {
        if (candidateList != null) {
            throw new IllegalStateException("Cannot add more alt contexts");
        }

        for (final AltContext altContext : altContexts) {
            if (altContext.primaryReadContext().readContext().isComplete()) {
                candidateMap.computeIfAbsent(altContext, x -> new Candidate(altContext)).update(altContext);
            }
        }
    }

    @NotNull
    public List<Candidate> candidates() {
        if (candidateList == null) {
            final VariantHotspotComparator variantHotspotComparator = new VariantHotspotComparator();
            candidateList = Lists.newArrayList(candidateMap.values());
            candidateList.sort((o1, o2) -> variantHotspotComparator.compare(o1.variant(), o2.variant()));
        }

        return candidateList;
    }

}
