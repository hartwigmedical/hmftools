package com.hartwig.hmftools.amber;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.jetbrains.annotations.NotNull;

public class AmberHetNormalEvidence
{
    private final Map<String, Collection<BaseDepth>> mMap;
    private final BaseDepthIntersectFilter mIntersectFilter;

    public AmberHetNormalEvidence()
    {
        mMap = Maps.newHashMap();
        mIntersectFilter = new BaseDepthIntersectFilter();
    }

    @NotNull
    public ListMultimap<Chromosome, AmberSite> intersection()
    {
        return mIntersectFilter.sites();
    }

    @NotNull
    public Set<String> samples()
    {
        return mMap.keySet();
    }

    public Collection<BaseDepth> evidence(String sample)
    {
        return mMap.get(sample);
    }

    public void add(@NotNull final String sample, @NotNull Collection<BaseDepth> baseDepths)
    {
        mMap.put(sample, baseDepths);
        mIntersectFilter.additional(baseDepths);
    }

    @NotNull
    public Predicate<BaseDepth> intersectionFilter()
    {
        return mIntersectFilter;
    }
}
