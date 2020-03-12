package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.context.RefContextFixedFactoryForSample;

import org.jetbrains.annotations.NotNull;

public class AltContextCollection {

    private final SageConfig sageConfig;
    private final Set<VariantHotspot> hotspots;
    private final ListMultimap<VariantHotspot, AltContext> map = ArrayListMultimap.create();
    private final Comparator<AltContext> comparator;

    public AltContextCollection(@NotNull final SageConfig config, @NotNull final String primarySample, @NotNull final List<VariantHotspot> hotspots) {
        this.sageConfig = config;
        this.hotspots = Sets.newHashSet(hotspots);
        this.comparator = (o1, o2) -> {
            if (o1.sample().equals(primarySample)) {
                return -1;
            }

            if (o2.sample().equals(primarySample)) {
                return 1;
            }

            return o1.sample().compareTo(o2.sample());
        };
    }

    public void addRefContexts(@NotNull final Collection<RefContext> refContexts) {
        for (RefContext refContext : refContexts) {
            addAltContexts(refContext.alts());
        }
    }

    private void addAltContexts(@NotNull final Collection<AltContext> altContexts) {
        for (final AltContext altContext : altContexts) {
            if (hotspots.contains(altContext)) {
                map.put(altContext, altContext);
            }
        }
    }

    @NotNull
    public List<AltContext> altContexts(@NotNull final VariantHotspot hotspot) {
        assert (map.containsKey(hotspot));
        final List<AltContext> result = map.get(hotspot);
        result.sort(comparator);
        return map.get(hotspot);
    }

    @NotNull
    public RefContextFixedFactoryForSample createFactory(@NotNull final Predicate<AltContext> anyPredicate) {
        return new RefContextFixedFactoryForSample(sageConfig, anyMatch(anyPredicate));
    }

    @NotNull
    private Map<VariantHotspot, Candidate> anyMatch(@NotNull final Predicate<AltContext> predicate) {
        final Map<VariantHotspot, Candidate> result = Maps.newHashMap();

        for (final VariantHotspot variant : map.keySet()) {
            final List<AltContext> list = map.get(variant);
            if (list.stream().anyMatch(predicate)) {
                Candidate candidate = new Candidate(variant);
                list.forEach(candidate::update);
                result.put(candidate.hotspot(), candidate);
            }

        }

        return result;
    }

}
