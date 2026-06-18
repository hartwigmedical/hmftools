package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConfig.TARS_LOGGER;
import static com.hartwig.hmftools.tars.liftback.ChunkProducer.END_OF_STREAM;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import com.hartwig.hmftools.tars.liftback.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.tars.liftback.rescue.RefSequenceSource;
import com.hartwig.hmftools.tars.liftback.rescue.RescueStatistics;
import com.hartwig.hmftools.tars.liftback.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.tars.liftback.tailextend.TailExtensionStatistics;
import com.hartwig.hmftools.tars.liftback.tailextend.TerminalMicroJunctionCollapser;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

// Drains read-name chunks from the shared queue, transforms each contiguous name-group with a per-group
// LiftedMateInfoCache, and writes lifted records to its own unsorted shard. Owns its engines and ref-source
// handle so nothing serialises across workers. Exits on END_OF_STREAM.
public class LiftBackWorker extends Thread
{
    private final BlockingQueue<List<SAMRecord>> mQueue;
    private final LiftBackStats mStats;
    private final LiftBackGroupProcessor mProcessor;
    private final SAMFileWriter mShardWriter;
    private final LiftBackWriter mTsvWriter; // nullable: headerless per-worker TSV shard, only when enabled

    private final JunctionRescueResolver mRescueResolver;
    private final SoftclipTailExtender mSoftclipExtender;
    private final TerminalMicroJunctionCollapser mTerminalCollapser;
    private final JunctionCanonicalizer mJunctionCanonicalizer;

    private final ExcludedRegions mExcludedRegions; // nullable: passed to the processor for post-lift exclusion

    public LiftBackWorker(
            final BlockingQueue<List<SAMRecord>> queue, final LiftBackResources resources,
            final SAMFileHeader header, final String shardBam, final String tsvAShard, final String tsvBShard)
    {
        mQueue = queue;
        mStats = new LiftBackStats();
        mExcludedRegions = resources.ExcludedRegions;

        try
        {
            mTsvWriter = tsvAShard != null ? new LiftBackWriter(tsvAShard, tsvBShard, false) : null;
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to open liftback TSV shard: " + e, e);
        }

        final RefSequenceSource refSource = resources.openRefSource();

        mRescueResolver = new JunctionRescueResolver(resources.JunctionIndex, refSource, resources.Rescue);
        mSoftclipExtender = new SoftclipTailExtender(refSource, resources.JunctionIndex, resources.TailExtension);
        mTerminalCollapser = refSource != null
                ? new TerminalMicroJunctionCollapser(refSource, resources.TerminalAnchor) : null;
        mJunctionCanonicalizer = refSource != null
                ? new JunctionCanonicalizer(refSource, JunctionCanonicalizer.DEFAULT_MAX_SHIFT) : null;

        mProcessor = new LiftBackGroupProcessor(
                resources.Resolver, mRescueResolver, mSoftclipExtender, mTerminalCollapser, mJunctionCanonicalizer,
                refSource, mExcludedRegions, mStats);

        mShardWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(shardBam));
    }

    public LiftBackStats liftBackStats() { return mStats; }
    public RescueStatistics rescueStatistics() { return mRescueResolver != null ? mRescueResolver.statistics() : null; }
    public TailExtensionStatistics tailExtStatistics() { return mSoftclipExtender != null ? mSoftclipExtender.statistics() : null; }
    public long collapsedLeading() { return mTerminalCollapser != null ? mTerminalCollapser.collapsedLeading() : 0; }
    public long collapsedTrailing() { return mTerminalCollapser != null ? mTerminalCollapser.collapsedTrailing() : 0; }
    public long junctionsCanonicalized() { return mJunctionCanonicalizer != null ? mJunctionCanonicalizer.junctionsShifted() : 0; }
    public long overCapUnmapped() { return mProcessor.overCapUnmapped(); }

    @Override
    public void run()
    {
        try
        {
            while(true)
            {
                final List<SAMRecord> chunk = mQueue.take();
                if(chunk == END_OF_STREAM)
                    break;
                processChunk(chunk);
            }
        }
        catch(InterruptedException e)
        {
            Thread.currentThread().interrupt();
            TARS_LOGGER.error("liftback worker interrupted: {}", e.toString());
            System.exit(1);
        }
        catch(Exception e)
        {
            TARS_LOGGER.error("liftback worker failed: {}", e.toString());
            System.exit(1);
        }
        finally
        {
            mShardWriter.close();
            if(mTsvWriter != null)
            {
                try
                {
                    mTsvWriter.close();
                }
                catch(IOException e)
                {
                    TARS_LOGGER.warn("failed to close liftback TSV shard: {}", e.toString());
                }
            }
        }
    }

    private void processChunk(final List<SAMRecord> chunk)
    {
        final List<SAMRecord> group = new ArrayList<>();
        String currentName = null;

        for(final SAMRecord record : chunk)
        {
            final String name = record.getReadName();
            if(currentName != null && !name.equals(currentName))
            {
                processNameGroup(group);
                group.clear();
            }
            group.add(record);
            currentName = name;
        }

        if(!group.isEmpty())
            processNameGroup(group);
    }

    private void processNameGroup(final List<SAMRecord> group)
    {
        // exclusion is applied post-lift inside the processor (on lifted genomic coords): a primary lifting into
        // an excluded region is unmapped REDUX-style, a supp dropped. Pre-lift can't test tx-contig reads, whose
        // input coords are chrN_tx, against the genomic region list.
        mProcessor.processNameGroup(group, new LiftedMateInfoCache(), this::write);
    }

    public long excludedReads() { return mProcessor.excludedReads(); }

    private void write(final SAMRecord record, final LiftBackResult result)
    {
        mShardWriter.addAlignment(record);
        if(mTsvWriter == null)
            return;
        try
        {
            mTsvWriter.write(record, result);
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to write liftback TSV shard: " + e, e);
        }
    }
}
