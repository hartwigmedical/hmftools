package com.hartwig.hmftools.amber;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthIntersectFilter;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.jetbrains.annotations.NotNull;

class AmberHetNormalEvidence {

    private final Map<String, Collection<BaseDepth>> map = Maps.newHashMap();
    private final BaseDepthIntersectFilter intersectFilter = new BaseDepthIntersectFilter();

    @NotNull
    public ListMultimap<Chromosome, AmberSite> intersection() {
        return intersectFilter.sites();
    }

    @NotNull
    public Set<String> samples() {
        return map.keySet();
    }

    public Collection<BaseDepth> evidence(String sample) {
        return map.get(sample);
    }

    public void add(@NotNull final String sample, @NotNull Collection<BaseDepth> baseDepths) {
        map.put(sample, baseDepths);
        intersectFilter.additional(baseDepths);
    }

    @NotNull
    public Predicate<BaseDepth> intersectionFilter() {
        return intersectFilter;
    }
}
