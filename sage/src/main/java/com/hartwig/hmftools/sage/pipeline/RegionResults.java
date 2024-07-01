package com.hartwig.hmftools.sage.pipeline;

import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.evidence.EvidenceStats;
import com.hartwig.hmftools.sage.sync.FragmentSyncType;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

public class RegionResults
{
    private final VcfWriter mVcfWriter;
    private int mTotalReads;
    private int mCandidates;
    private int mTotaVariants;
    private final List<PerformanceCounter> mPerfCounters;
    private final int[] mSyncCounts;
    private final EvidenceStats mEvidenceStats;

    public RegionResults(final VcfWriter vcfWriter)
    {
        mVcfWriter = vcfWriter;
        mCandidates = 0;
        mTotalReads = 0;
        mTotaVariants = 0;
        mPerfCounters = Lists.newArrayList();
        mSyncCounts = new int[FragmentSyncType.values().length];
        mEvidenceStats = new EvidenceStats();
    }

    public synchronized void addCandidates(int candidateCount)
    {
        mCandidates += candidateCount;
    }

    public synchronized void addFinalVariants(final int taskId, final List<SageVariant> variants)
    {
        mTotaVariants += variants.size();

        if(mVcfWriter != null)
            mVcfWriter.writeVariants(taskId, variants);
    }

    public synchronized void addTotalReads(int totalReads)
    {
        mTotalReads += totalReads;
    }

    public synchronized void addPerfCounters(final List<PerformanceCounter> perfCounters)
    {
        if(mPerfCounters.isEmpty())
        {
            mPerfCounters.addAll(perfCounters);
        }
        else
        {
            for(int j = 0; j < perfCounters.size(); ++j)
            {
                mPerfCounters.get(j).merge(perfCounters.get(j));
            }
        }
    }

    public synchronized void addSynCounts(final int[] counts)
    {
        Arrays.stream(FragmentSyncType.values()).forEach(x -> mSyncCounts[x.ordinal()] += counts[x.ordinal()]);
    }

    public synchronized void addEvidenceStats(final EvidenceStats stats)
    {
        mEvidenceStats.merge(stats);
    }

    public int totalReads() { return mTotalReads; }
    public int totalCandidates() { return mCandidates; }
    public int totalVariants() { return mTotaVariants; }

    public void logPerfCounters()
    {
        mPerfCounters.forEach(x -> x.logStats());
    }
    public EvidenceStats evidenceStats() { return mEvidenceStats; }

    public void logSynCounts()
    {
        if(!SG_LOGGER.isDebugEnabled())
            return;

        StringJoiner sj = new StringJoiner(", ");
        Arrays.stream(FragmentSyncType.values())
                .filter(x -> mSyncCounts[x.ordinal()] > 0)
                .forEach(x -> sj.add(format("%s=%d", x, mSyncCounts[x.ordinal()])));

        if(sj.length() > 0)
        {
            SG_LOGGER.debug("fragment sync counts: {}", sj);
        }
    }
}
