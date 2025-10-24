package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;
import com.hartwig.hmftools.sage.candidate.RefContextCache;
import com.hartwig.hmftools.sage.common.SamSlicerInterface;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.RefSequence;

import htsjdk.samtools.SAMRecord;

public class CandidateEvidence
{
    private final SageConfig mConfig;
    private final List<SimpleVariant> mHotspots;
    private final List<BaseRegion> mPanel;

    private int mTotalReadsProcessed;

    public CandidateEvidence(
            final SageConfig config, final List<SimpleVariant> hotspots, final List<BaseRegion> panel)
    {
        mConfig = config;
        mPanel = panel;
        mHotspots = hotspots;

        mTotalReadsProcessed = 0;
    }

    public int totalReadsProcessed() { return mTotalReadsProcessed; }

    public List<AltContext> readBam(final SamSlicerInterface samSlicer, final RefSequence refSequence, final ChrBaseRegion bounds)
    {
        final RefContextCache refContextCache = new RefContextCache(mConfig, mHotspots, mPanel);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(mConfig, bounds, refSequence, refContextCache, mHotspots);

        final Consumer<SAMRecord> consumer = record ->
        {
            try
            {
                refContextConsumer.processRead(record);
            }
            catch(Exception e)
            {
                SG_LOGGER.error("error processing read({} {}:{}-{}) error: {}",
                        record.getReadName(), record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd(), e.toString());
                e.printStackTrace();
                System.exit(1);
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
