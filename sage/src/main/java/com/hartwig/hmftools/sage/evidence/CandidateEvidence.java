package com.hartwig.hmftools.sage.evidence;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;
import com.hartwig.hmftools.sage.candidate.RefContextCache;
import com.hartwig.hmftools.sage.common.SamSlicerInterface;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.coverage.GeneCoverage;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.SAMRecord;

public class CandidateEvidence
{
    private final SageConfig mConfig;
    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanel;
    private final Coverage mCoverage;

    private int mTotalReadsProcessed;

    public CandidateEvidence(
            final SageConfig config, final List<VariantHotspot> hotspots, final List<BaseRegion> panel, final Coverage coverage)
    {
        mConfig = config;
        mPanel = panel;
        mHotspots = hotspots;
        mCoverage = coverage;

        mTotalReadsProcessed = 0;
    }

    public int totalReadsProcessed() { return mTotalReadsProcessed; }

    public List<AltContext> readBam(
            final String sample, final SamSlicerInterface samSlicer, final RefSequence refSequence, final ChrBaseRegion bounds)
    {
        final List<GeneCoverage> geneCoverage = mCoverage.coverage(sample, bounds.Chromosome);
        final RefContextCache refContextCache = new RefContextCache(mConfig, mHotspots, mPanel);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(mConfig, bounds, refSequence, refContextCache, mHotspots);

        final Consumer<SAMRecord> consumer = record ->
        {
            refContextConsumer.processRead(record);

            if(!geneCoverage.isEmpty())
            {
                geneCoverage.forEach(x -> x.processRead(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd()));
            }
        };

        List<AltContext> altContexts = readBam(samSlicer, consumer, refContextCache);

        mTotalReadsProcessed += refContextConsumer.getReadCount();

        return altContexts;
    }

    private List<AltContext> readBam(
            final SamSlicerInterface samSlicer, final Consumer<SAMRecord> recordConsumer, final RefContextCache refContextCache)
    {
        final List<AltContext> altContexts = Lists.newArrayList();

        samSlicer.slice(recordConsumer);

        // add all valid alt contexts
        altContexts.addAll(refContextCache.altContexts());

        return altContexts;
    }
}
