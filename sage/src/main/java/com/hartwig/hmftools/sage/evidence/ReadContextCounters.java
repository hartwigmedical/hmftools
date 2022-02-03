package com.hartwig.hmftools.sage.evidence;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.FilterConfig;

public class ReadContextCounters
{
    private final List<Candidate> mCandidates;

    private final Map<VariantHotspot,List<ReadContextCounter>> mVariantReadCounters;

    public ReadContextCounters(final String primarySample, final List<Candidate> candidates)
    {
        mCandidates = candidates;
        mVariantReadCounters = Maps.newHashMap();
    }

    public int variantCount() { return mVariantReadCounters.size(); }

    public List<ReadContextCounter> getVariantReadCounters(final VariantHotspot variant)
    {
        return mVariantReadCounters.get(variant);
    }

    public void addCounters(final List<ReadContextCounter> readCounters, int expectedCount)
    {
        for(ReadContextCounter readCounter : readCounters)
        {
            List<ReadContextCounter> variantReadCounters = mVariantReadCounters.get(readCounter.variant());
            if(variantReadCounters == null)
            {
                variantReadCounters = Lists.newArrayListWithExpectedSize(expectedCount);
                mVariantReadCounters.put(readCounter.variant(), variantReadCounters);

            }

            variantReadCounters.add(readCounter);
        }
    }

    public List<Candidate> filterCandidates(final FilterConfig filterConfig)
    {
        final List<Candidate> result = Lists.newArrayList();

        for(Candidate candidate : mCandidates)
        {
            List<ReadContextCounter> readCounters = mVariantReadCounters.get(candidate.variant());

            if(readCounters != null && readCounters.stream().anyMatch(x -> filterConfig.passesHardFilters(x)))
                result.add(candidate);
        }
        return result;
    }
}
