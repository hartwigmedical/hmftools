package com.hartwig.hmftools.sage.read;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;

import org.jetbrains.annotations.NotNull;

public class ReadContextCounters
{
    private final Comparator<ReadContextCounter> mComparator;
    private final List<Candidate> mCandidates;
    private final ListMultimap<VariantHotspot, ReadContextCounter> mMap = ArrayListMultimap.create();

    public ReadContextCounters(@NotNull final String primarySample, @NotNull final List<Candidate> candidates)
    {
        mCandidates = candidates;
        mComparator = (o1, o2) ->
        {
            if(o1.Sample.equals(primarySample))
                return -1;

            if(o2.Sample.equals(primarySample))
                return 1;

            return o1.Sample.compareTo(o2.Sample);
        };
    }

    @NotNull
    public List<ReadContextCounter> readContextCounters(@NotNull final VariantHotspot variant)
    {
        assert (mMap.containsKey(variant));
        final List<ReadContextCounter> result = mMap.get(variant);
        result.sort(mComparator);
        return result;
    }

    public void addCounters(Collection<ReadContextCounter> counters)
    {
        for(ReadContextCounter counter : counters)
        {
            mMap.put(counter.Variant, counter);
        }
    }

    @NotNull
    public List<Candidate> candidates(@NotNull final Predicate<ReadContextCounter> anyPredicate)
    {
        final List<Candidate> result = Lists.newArrayList();
        for(Candidate candidate : mCandidates)
        {
            List<ReadContextCounter> counters = mMap.get(candidate.variant());
            if(counters != null && counters.stream().anyMatch(anyPredicate))
            {
                result.add(candidate);
            }
        }
        return result;
    }
}
