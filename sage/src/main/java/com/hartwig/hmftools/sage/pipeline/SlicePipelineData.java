package com.hartwig.hmftools.sage.pipeline;

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
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.NormalRefContextCandidates;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefContextCandidates;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.jetbrains.annotations.NotNull;

public class SlicePipelineData {

    private final String normalSample;

    private final Set<VariantHotspot> allHotspots = Sets.newHashSet();
    private final List<Map<VariantHotspot, AltContext>> altContextMap = new ArrayList<>();
    private final Map<Long, RefContext> normalMap = Maps.newHashMap();
    private final SageVariantFactory variantFactory;

    SlicePipelineData(final String normalSample, final int tumorSampleSize, final SageVariantFactory variantFactory) {
        this.normalSample = normalSample;
        this.variantFactory = variantFactory;
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
                readContextCounters.sort(Comparator.comparingInt(ReadContextCounter::support).reversed());
                altContext.setPrimaryReadContext(new ReadContextCounter(hotspot, readContextCounters.get(0).readContext()));
            }
        }

        return candidates;
    }

    public void addNormal(@NotNull final List<RefContext> altContexts) {
        altContexts.forEach(x -> normalMap.put(x.position(), x));
    }

    @NotNull
    public List<SageVariant> results() {
        List<SageVariant> result = Lists.newArrayList();

        final Comparator<VariantHotspot> hotspotComparator = (o1, o2) -> {
            int standardCompare = o1.compareTo(o2);
            if (standardCompare != 0) {
                return standardCompare;
            }

            int o1Length = Math.max(o1.ref().length(), o1.alt().length());
            int o2Length = Math.max(o2.ref().length(), o2.alt().length());
            int lengthCompare = Integer.compare(o1Length, o2Length);
            if (lengthCompare != 0) {
                return lengthCompare;
            }

            int refCompare = o1.ref().compareTo(o2.ref());
            if (refCompare != 0) {
                return refCompare;
            }

            return o1.alt().compareTo(o2.alt());
        };

        final List<VariantHotspot> sortedHotspots = allHotspots.stream().sorted(hotspotComparator).collect(Collectors.toList());

        for (VariantHotspot sortedHotspot : sortedHotspots) {

            final RefContext normalRefContext = normalMap.get(sortedHotspot.position());
            final AltContext normalAltContext = normalRefContext == null
                    ? new AltContext(normalSample, sortedHotspot)
                    : normalRefContext.altContext(sortedHotspot.ref(), sortedHotspot.alt());

            final List<AltContext> tumorAltContexts = new ArrayList<>(altContextMap.size() + 1);
            for (Map<VariantHotspot, AltContext> variantHotspotAltContextMap : altContextMap) {
                AltContext tumorContext = variantHotspotAltContextMap.get(sortedHotspot);
                if (tumorContext != null) {
                    tumorAltContexts.add(tumorContext);
                }
            }

            result.add(variantFactory.create(normalAltContext, tumorAltContexts));
        }

        return result;
    }

}
