package com.hartwig.hmftools.sage.pipeline;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.vcf.VcfWriter;

public class RegionResults
{
    private final VcfWriter mVcfWriter;
    private int mTotalReads;
    private int mTotaVariants;
    private final List<PerformanceCounter> mPerfCounters;

    public RegionResults(final VcfWriter vcfWriter)
    {
        mVcfWriter = vcfWriter;
        mTotalReads = 0;
        mTotaVariants = 0;
        mPerfCounters = Lists.newArrayList();
    }

    public synchronized void addFinalVariants(final int taskId, final List<SageVariant> variants)
    {
        mTotaVariants += variants.size();
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

    public int totalReads() { return mTotalReads; }
    public int totalVariants() { return mTotaVariants; }

    public void logPerfCounters()
    {
        if(SG_LOGGER.isDebugEnabled())
        {
            mPerfCounters.forEach(x -> x.logStats());
        }
    }

}
