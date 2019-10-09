package com.hartwig.hmftools.sage;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.NormalRefContextCandidates;
import com.hartwig.hmftools.sage.context.ReadContextCounter;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefContextCandidates;

import org.jetbrains.annotations.NotNull;

public class SagePipelineData {

    private final String normalSample;

    private final Set<VariantHotspot> allHotspots = Sets.newHashSet();
    private final List<Map<VariantHotspot, AltContext>> altContextMap = new ArrayList<>();
    private final Map<Long, RefContext> normalMap = Maps.newHashMap();

    SagePipelineData(final String normalSample, final int tumorSampleSize) {
        this.normalSample = normalSample;
        for (int i = 0; i < tumorSampleSize; i++) {
            altContextMap.add(new HashMap<>());
        }
    }

    public void addTumor(int tumorSample, @NotNull final List<AltContext> altContexts) {
        final Map<VariantHotspot, AltContext> map = altContextMap.get(tumorSample);
        altContexts.forEach(x -> {
            allHotspots.add(x);
            map.put(x, x);
        });

    }

    @NotNull
    public RefContextCandidates normalCandidates() {

        NormalRefContextCandidates candidates = new NormalRefContextCandidates(normalSample);

        List<VariantHotspot> sortedHotspots =
                allHotspots.stream().sorted(Comparator.comparingLong(GenomePosition::position)).collect(Collectors.toList());

        for (VariantHotspot hotspot : sortedHotspots) {
            final RefContext refContext = candidates.add(hotspot.chromosome(), hotspot.position());
            final AltContext altContext = refContext.altContext(hotspot.ref(), hotspot.alt());

            final List<ReadContextCounter> readContextCounters = Lists.newArrayList();
            for (Map<VariantHotspot, AltContext> variantHotspotAltContextMap : altContextMap) {
                AltContext tumorContext = variantHotspotAltContextMap.get(hotspot);
                if (tumorContext != null) {
                    readContextCounters.add(tumorContext.primaryReadContext());
                }
            }

            if (!readContextCounters.isEmpty()) {
                readContextCounters.sort(Comparator.comparingInt(ReadContextCounter::full).reversed());
                altContext.setPrimaryReadContext(new ReadContextCounter(hotspot, readContextCounters.get(0).readContext()));
            }
        }

        return candidates;
    }

    public void addNormal(@NotNull final List<RefContext> altContexts) {
        altContexts.forEach(x -> normalMap.put(x.position(), x));
    }

    @NotNull
    public List<List<AltContext>> altContexts() {
        List<List<AltContext>> result = Lists.newArrayList();

        final Comparator<VariantHotspot> hotspotComparator = (o1, o2) -> {
            int standardCompare = o1.compareTo(o2);
            if (standardCompare != 0) {
                return standardCompare;
            }

            int refCompare = o1.ref().compareTo(o2.ref());
            if (refCompare != 0) {
                return refCompare;
            }

            return o1.alt().compareTo(o2.alt());
        };

        final List<VariantHotspot> sortedHotspots = allHotspots.stream().sorted(hotspotComparator).collect(Collectors.toList());

        for (VariantHotspot sortedHotspot : sortedHotspots) {

            final List<AltContext> altContexts = new ArrayList<>(altContextMap.size() + 1);
            final RefContext normalRefContext = normalMap.get(sortedHotspot.position());
            if (normalRefContext != null) {
                altContexts.add(normalRefContext.altContext(sortedHotspot.ref(), sortedHotspot.alt()));
            } else {
                altContexts.add(new AltContext(normalSample, sortedHotspot));
            }

            for (Map<VariantHotspot, AltContext> variantHotspotAltContextMap : altContextMap) {
                AltContext tumorContext = variantHotspotAltContextMap.get(sortedHotspot);
                if (tumorContext != null) {
                    altContexts.add(tumorContext);
                }
            }

            result.add(altContexts);
        }

        return result;
    }

}
