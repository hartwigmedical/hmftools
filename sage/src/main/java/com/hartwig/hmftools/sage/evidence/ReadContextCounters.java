package com.hartwig.hmftools.sage.evidence;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.FilterConfig;

public class ReadContextCounters
{
    private final List<Candidate> mCandidates;

    // multiple read counters exist to support multiple samples
    private final List<List<ReadContextCounter>> mSampleCandidateReadCounters;
    private final List<Integer> mFilteredCandidateIndex; // index of the filtered candidates into the original/full list

    // private final Map<VariantHotspot,List<ReadContextCounter>> mVariantReadCounters;

    public ReadContextCounters(final List<Candidate> candidates)
    {
        mCandidates = candidates;
        mSampleCandidateReadCounters = Lists.newArrayList();
        // mVariantReadCounters = Maps.newHashMap();
        mFilteredCandidateIndex = Lists.newArrayList();
    }

    public int candidateCount() { return mCandidates.size(); }
    // public int variantCount() { return mVariantReadCounters.size(); }

    public List<ReadContextCounter> getReadCounters(final int candidateIndex)
    {
        if(candidateIndex < 0 ||  candidateIndex >= mCandidates.size())
            return null;

        return mSampleCandidateReadCounters.get(candidateIndex);
    }

    public List<ReadContextCounter> getFilteredReadCounters(final int filteredCandidateIndex)
    {
        if(filteredCandidateIndex < 0 ||  filteredCandidateIndex >= mFilteredCandidateIndex.size())
            return null;

        int candidateIndex = mFilteredCandidateIndex.get(filteredCandidateIndex);
        return mSampleCandidateReadCounters.get(candidateIndex);
    }

    public List<ReadContextCounter> getReadCounters(final Candidate candidate)
    {
        for(int i = 0; i < mCandidates.size(); ++i)
        {
            if(mCandidates.get(i) == candidate)
                return mSampleCandidateReadCounters.get(i);
        }

        return null;
    }

    public List<ReadContextCounter> getVariantReadCounters(final VariantHotspot variant)
    {
        final VariantHotspotComparator variantComparator = new VariantHotspotComparator();

        List<ReadContextCounter> readCounters = Lists.newArrayList();

        mSampleCandidateReadCounters.stream()
                .filter(x -> variantComparator.compare(x.get(0).variant(), variant) == 0)
                .forEach(x -> readCounters.addAll(x));

        return readCounters;
    }

    public void addCounters(final List<ReadContextCounter> sampleReadCounters, int expectedCount)
    {
        if(sampleReadCounters.size() != mCandidates.size())
        {
            throw new IllegalStateException("non-matching read-counters with candidates");
        }

        boolean firstSample = mSampleCandidateReadCounters.isEmpty();

        for(int i = 0; i < mCandidates.size(); ++i)
        {
            List<ReadContextCounter> readCounters;

            if(firstSample)
            {
                readCounters = Lists.newArrayListWithExpectedSize(expectedCount);
                mSampleCandidateReadCounters.add(readCounters);
            }
            else
            {
                readCounters = mSampleCandidateReadCounters.get(i);
            }

            readCounters.add(sampleReadCounters.get(i));
        }

        /*
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
        */
    }

    public List<Candidate> filterCandidates(final FilterConfig filterConfig)
    {
        final List<Candidate> validCandidates = Lists.newArrayList();

        for(int i = 0; i < mCandidates.size(); ++i)
        {
            List<ReadContextCounter> readCounters = mSampleCandidateReadCounters.get(i);

            if(readCounters.stream().anyMatch(x -> filterConfig.passesHardFilters(x)))
            {
                Candidate candidate = mCandidates.get(i);

                validCandidates.add(mCandidates.get(i));
                mFilteredCandidateIndex.add(i);
            }
        }

        /*
        for(Candidate candidate : mCandidates)
        {
            List<ReadContextCounter> readCounters = mVariantReadCounters.get(candidate.variant());

            if(readCounters != null && readCounters.stream().anyMatch(x -> filterConfig.passesHardFilters(x)))
                validCandidates.add(candidate);
        }
        */

        return  validCandidates;
    }
}
