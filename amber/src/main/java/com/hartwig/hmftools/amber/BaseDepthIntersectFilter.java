package com.hartwig.hmftools.amber;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;

public class BaseDepthIntersectFilter implements Predicate<BaseDepth>
{
    private boolean mAdditional;
    private final Set<AmberSite> mIntersection;

    public BaseDepthIntersectFilter()
    {
        mAdditional = false;
        mIntersection = Sets.newHashSet();
    }

    @NotNull
    public ListMultimap<Chromosome, AmberSite> sites()
    {
        ListMultimap<Chromosome, AmberSite> result = ArrayListMultimap.create();

        for(AmberSite amberSite : mIntersection)
        {
            result.put(HumanChromosome.fromString(amberSite.chromosome()), amberSite);
        }

        for(Chromosome chromosome : result.keySet())
        {
            List<AmberSite> list = result.get(chromosome);
            Collections.sort(list);
        }

        return result;
    }

    public void additional(@NotNull final Collection<BaseDepth> depth)
    {
        final Set<AmberSite> set = asSet(depth);
        if(mAdditional)
        {
            mIntersection.retainAll(set);
        }
        else
        {
            mAdditional = true;
            mIntersection.addAll(set);
        }
    }

    public int size()
    {
        return mIntersection.size();
    }

    @Override
    public boolean test(final BaseDepth baseDepth)
    {
        return !mAdditional || mIntersection.contains(AmberSiteFactory.asSite(baseDepth));
    }

    @NotNull
    private static Set<AmberSite> asSet(@NotNull final Collection<BaseDepth> depth)
    {
        return depth.stream().map(AmberSiteFactory::asSite).collect(Collectors.toSet());
    }

}
