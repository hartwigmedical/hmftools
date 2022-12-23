package com.hartwig.hmftools.sage.pipeline;

import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.evidence.SyncFragmentType;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

public class RegionResults
{
    private final VcfWriter mVcfWriter;
    private int mTotalReads;
    private int mTotaVariants;
    private int mMaxMemoryUsage;
    private final List<PerformanceCounter> mPerfCounters;
    private final int[] mSyncCounts;

    public RegionResults(final VcfWriter vcfWriter)
    {
        mVcfWriter = vcfWriter;
        mTotalReads = 0;
        mTotaVariants = 0;
        mMaxMemoryUsage = 0;
        mPerfCounters = Lists.newArrayList();
        mSyncCounts = new int[SyncFragmentType.values().length];
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

    public synchronized void addMaxMemory(int maxMemory)
    {
        mMaxMemoryUsage = max(mMaxMemoryUsage, maxMemory);
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
        Arrays.stream(SyncFragmentType.values()).forEach(x -> mSyncCounts[x.ordinal()] += counts[x.ordinal()]);
    }

    public int totalReads() { return mTotalReads; }
    public int totalVariants() { return mTotaVariants; }
    public int maxMemoryUsage() { return mMaxMemoryUsage; }

    public void logPerfCounters()
    {
        mPerfCounters.forEach(x -> x.logStats());
    }

    public void logSynCounts()
    {
        StringJoiner sj = new StringJoiner(", ");
        Arrays.stream(SyncFragmentType.values()).forEach(x -> sj.add(format("%s=%d", x, mSyncCounts[x.ordinal()])));
        SG_LOGGER.info("fragment sync counts: {}", sj);
    }
}
