package com.hartwig.hmftools.sage.pipeline;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.NormalRefContextCandidates;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.jetbrains.annotations.NotNull;

class SomaticPipelineData {

    private final SageConfig config;
    private final Set<VariantHotspot> allHotspots = Sets.newHashSet();
    private final List<Map<Long, RefContext>> normalMap = new ArrayList<>();
    private final List<Map<VariantHotspot, AltContext>> tumorMap = new ArrayList<>();
    private final SageVariantFactory variantFactory;

    SomaticPipelineData(final SageConfig config, final SageVariantFactory variantFactory) {
        this.config = config;
        this.variantFactory = variantFactory;
        for (int i = 0; i < config.tumor().size(); i++) {
            tumorMap.add(new HashMap<>());
        }
        for (int i = 0; i < config.reference().size(); i++) {
            normalMap.add(new HashMap<>());
        }
    }

    public void addTumor(int sampleNumber, @NotNull final List<AltContext> altContexts) {
        final Map<VariantHotspot, AltContext> map = tumorMap.get(sampleNumber);
        altContexts.forEach(x -> {
            allHotspots.add(x);
            map.put(x, x);
        });
    }

    public void addNormal(int sampleNumber, @NotNull final List<RefContext> refContexts) {
        final Map<Long, RefContext> map = normalMap.get(sampleNumber);
        refContexts.forEach(x -> map.put(x.position(), x));
    }

    @NotNull
    public NormalRefContextCandidates normalCandidates(final String sample) {

        NormalRefContextCandidates candidates = new NormalRefContextCandidates(sample);

        List<VariantHotspot> sortedHotspots =
                allHotspots.stream().sorted(Comparator.comparingLong(GenomePosition::position)).collect(Collectors.toList());

        for (VariantHotspot hotspot : sortedHotspots) {
            final RefContext refContext = candidates.add(hotspot.chromosome(), hotspot.position());
            final AltContext altContext = refContext.altContext(hotspot.ref(), hotspot.alt());

            final List<ReadContextCounter> readContextCounters = Lists.newArrayList();
            for (Map<VariantHotspot, AltContext> variantHotspotAltContextMap : tumorMap) {
                AltContext tumorContext = variantHotspotAltContextMap.get(hotspot);
                if (tumorContext != null) {
                    readContextCounters.add(tumorContext.primaryReadContext());
                }
            }

            if (!readContextCounters.isEmpty()) {
                readContextCounters.sort(Comparator.comparingInt(ReadContextCounter::altSupport).reversed());
                altContext.setPrimaryReadContext(new ReadContextCounter(hotspot, readContextCounters.get(0).readContext()));
            }
        }

        return candidates;
    }

    @NotNull
    public List<SageVariant> results() {
        return results(variantFactory);
    }

    @NotNull
    public List<SageVariant> results(@NotNull final SageVariantFactory variantFactory) {
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

            final List<AltContext> normalAltContexts = new ArrayList<>(normalMap.size() + 1);
            for (int i = 0; i < normalMap.size(); i++) {
                final Map<Long, RefContext> contextMap = normalMap.get(i);
                final Optional<RefContext> normalRefContext = Optional.ofNullable(contextMap.get(sortedHotspot.position()));
                final Optional<AltContext> normalAltContext =
                        normalRefContext.map(x -> x.altContext(sortedHotspot.ref(), sortedHotspot.alt()));

                if (normalAltContext.isPresent()) {
                    normalAltContexts.add(normalAltContext.get());
                } else if (i == 0) {
                    // Primary normal cannot be absent
                    normalAltContexts.add(new AltContext(config.reference().get(i), sortedHotspot));
                }
            }

            final List<AltContext> tumorAltContexts = new ArrayList<>(tumorMap.size() + 1);
            for (Map<VariantHotspot, AltContext> variantHotspotAltContextMap : tumorMap) {
                AltContext tumorContext = variantHotspotAltContextMap.get(sortedHotspot);
                if (tumorContext != null) {
                    tumorAltContexts.add(tumorContext);
                }
            }

            result.add(variantFactory.create(normalAltContexts, tumorAltContexts));
        }

        return result;
    }

}
