package com.hartwig.hmftools.amber;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.jetbrains.annotations.NotNull;

public class HetNormalEvidence
{
    private final Map<String, Collection<PositionEvidence>> mMap;
    private final BaseDepthIntersectFilter mIntersectFilter;

    public HetNormalEvidence()
    {
        mMap = Maps.newHashMap();
        mIntersectFilter = new BaseDepthIntersectFilter();
    }

    public ListMultimap<Chromosome, AmberSite> intersection()
    {
        return mIntersectFilter.sites();
    }

    public Set<String> samples()
    {
        return mMap.keySet();
    }

    public Collection<PositionEvidence> evidence(String sample)
    {
        return mMap.get(sample);
    }

    public void add(@NotNull final String sample, @NotNull Collection<PositionEvidence> positionEvidences)
    {
        mMap.put(sample, positionEvidences);
        mIntersectFilter.additional(positionEvidences);
    }

    @NotNull
    public Predicate<PositionEvidence> intersectionFilter()
    {
        return mIntersectFilter;
    }
}
