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
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.context.FixedRefContextCandidatesFactory;
import com.hartwig.hmftools.sage.context.RefContext;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class AltContextCandidates {

    private final Set<VariantHotspot> hotspots;
    private final ListMultimap<VariantHotspot, AltContext> map = ArrayListMultimap.create();
    private final Comparator<AltContext> comparator;

    public AltContextCandidates(String primarySample, List<VariantHotspot> hotspots) {
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

    public void addRefCandidate(@NotNull final Collection<RefContext> refContexts) {
        for (RefContext refContext : refContexts) {
            addAltCandidate(refContext.alts());
        }
    }

    private void addAltCandidate(@NotNull final Collection<AltContext> altContexts) {
        for (final AltContext altContext : altContexts) {
            final VariantHotspot variant = ImmutableVariantHotspotImpl.builder().from(altContext).build();
            if (hotspots.contains(variant)) {
                map.put(variant, altContext);
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
    public FixedRefContextCandidatesFactory createFactory(@NotNull final Predicate<AltContext> anyPredicate) {
        return new FixedRefContextCandidatesFactory(anyMatch(anyPredicate));
    }

    @NotNull
    private Map<VariantHotspot, ReadContext> anyMatch(@NotNull final Predicate<AltContext> predicate) {
        final Map<VariantHotspot, ReadContext> result = Maps.newHashMap();

        for (final VariantHotspot variant : map.keySet()) {
            final List<AltContext> list = map.get(variant);
            if (list.stream().anyMatch(predicate)) {
                result.put(variant, list.get(0).primaryReadContext().readContext());
            }

        }

        return result;
    }

}
