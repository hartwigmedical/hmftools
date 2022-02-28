package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.CompletionException;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.candidate.AltContext;
import com.hartwig.hmftools.sage.candidate.RefContextConsumer;
import com.hartwig.hmftools.sage.candidate.RefContextCache;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.coverage.GeneCoverage;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.SamSlicer;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class CandidateEvidence
{
    private final SageConfig mConfig;
    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanel;
    private final ReferenceSequenceFile mRefGenome;
    private final Coverage mCoverage;

    private int mTotalReadsProcessed;

    public CandidateEvidence(
            final SageConfig config, final List<VariantHotspot> hotspots, final List<BaseRegion> panel,
            final ReferenceSequenceFile refGenome, final Coverage coverage)
    {
        mConfig = config;
        mPanel = panel;
        mHotspots = hotspots;
        mRefGenome = refGenome;
        mCoverage = coverage;

        mTotalReadsProcessed = 0;
    }

    public int totalReadsProcessed() { return mTotalReadsProcessed; }

    public List<AltContext> readBam(
            final String sample, final SamReader bamReader, final RefSequence refSequence, final ChrBaseRegion bounds)
    {
        final List<GeneCoverage> geneCoverage = mCoverage.coverage(sample, bounds.Chromosome);
        final RefContextCache refContextCache = new RefContextCache(mConfig, mHotspots, mPanel);
        final RefContextConsumer refContextConsumer = new RefContextConsumer(mConfig, bounds, refSequence, refContextCache, mHotspots);

        final Consumer<SAMRecord> consumer = record ->
        {
            refContextConsumer.accept(record);

            if(!geneCoverage.isEmpty())
            {
                geneCoverage.forEach(x -> x.processRead(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd()));
            }
        };

        List<AltContext> altContexts = readBam(bamReader, bounds, consumer, refContextCache);

        mTotalReadsProcessed += refContextConsumer.getReadCount();

        return altContexts;
    }

    @NotNull
    private List<AltContext> readBam(
            final SamReader bamReader, final ChrBaseRegion bounds, final Consumer<SAMRecord> recordConsumer, final RefContextCache refContextCache)
    {
        final List<AltContext> altContexts = Lists.newArrayList();

        final SamSlicer slicer = mConfig.PanelOnly ?
                new SamSlicer(0, bounds, mPanel) : new SamSlicer(0, bounds);

        slicer.slice(bamReader, recordConsumer);

        // add all valid alt contexts
        altContexts.addAll(refContextCache.altContexts());

        return altContexts;
    }
}
