package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.splice.ChunkProducer.END_OF_STREAM;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.BlockingQueue;

import com.hartwig.hmftools.redux.splice.rescue.JunctionRescueResolver;
import com.hartwig.hmftools.redux.splice.rescue.RefSequenceSource;
import com.hartwig.hmftools.redux.splice.rescue.RescueConfig;
import com.hartwig.hmftools.redux.splice.rescue.RescueStatistics;
import com.hartwig.hmftools.redux.splice.tailextend.SoftclipTailExtender;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionConfig;
import com.hartwig.hmftools.redux.splice.tailextend.TailExtensionStatistics;
import com.hartwig.hmftools.redux.splice.tailextend.TerminalMicroJunctionCollapser;

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
    private final LiftBackWriter mTsvWriter; // headerless per-worker records + alignments shard

    private final JunctionRescueResolver mRescueResolver;
    private final SoftclipTailExtender mSoftclipExtender;
    private final TerminalMicroJunctionCollapser mTerminalCollapser;
    private final JunctionCanonicalizer mJunctionCanonicalizer;

    private final ExcludedRegions mExcludedRegions; // nullable: drop fragments here before lifting
    private long mExcludedReads;

    public LiftBackWorker(
            final BlockingQueue<List<SAMRecord>> queue, final LiftBackResources resources,
            final SAMFileHeader header, final String shardBam, final String tsvAShard, final String tsvBShard)
    {
        mQueue = queue;
        mStats = new LiftBackStats();
        mExcludedRegions = resources.ExcludedRegions;

        try
        {
            mTsvWriter = new LiftBackWriter(tsvAShard, tsvBShard, false);
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to open liftback TSV shard: " + e, e);
        }

        final RefSequenceSource refSource =
                resources.RescueViaSupp || resources.ExtendSoftclipTails ? resources.openRefSource() : null;

        mRescueResolver = resources.RescueViaSupp
                ? new JunctionRescueResolver(resources.JunctionIndex, refSource, RescueConfig.enabledDefaults()) : null;
        mSoftclipExtender = resources.ExtendSoftclipTails
                ? new SoftclipTailExtender(refSource, resources.JunctionIndex, TailExtensionConfig.enabledDefaults()) : null;
        mTerminalCollapser = refSource != null
                ? new TerminalMicroJunctionCollapser(refSource, resources.TerminalAnchor) : null;
        mJunctionCanonicalizer = refSource != null
                ? new JunctionCanonicalizer(refSource, JunctionCanonicalizer.DEFAULT_MAX_SHIFT) : null;

        mProcessor = new LiftBackGroupProcessor(
                resources.Resolver, mRescueResolver, mSoftclipExtender, mTerminalCollapser, mJunctionCanonicalizer,
                refSource, resources.UnmapAboveNh, resources.UnmapBelowMapq, mStats);

        mShardWriter = new SAMFileWriterFactory().makeBAMWriter(header, false, new File(shardBam));
    }

    public LiftBackStats liftBackStats() { return mStats; }
    public RescueStatistics rescueStatistics() { return mRescueResolver != null ? mRescueResolver.statistics() : null; }
    public TailExtensionStatistics tailExtStatistics() { return mSoftclipExtender != null ? mSoftclipExtender.statistics() : null; }
    public long collapsedLeading() { return mTerminalCollapser != null ? mTerminalCollapser.collapsedLeading() : 0; }
    public long collapsedTrailing() { return mTerminalCollapser != null ? mTerminalCollapser.collapsedTrailing() : 0; }

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
            RD_LOGGER.error("liftback worker interrupted: {}", e.toString());
            System.exit(1);
        }
        catch(Exception e)
        {
            RD_LOGGER.error("liftback worker failed: {}", e.toString());
            System.exit(1);
        }
        finally
        {
            mShardWriter.close();
            try
            {
                mTsvWriter.close();
            }
            catch(IOException e)
            {
                RD_LOGGER.warn("failed to close liftback TSV shard: {}", e.toString());
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
        // pre-liftback filter: drop the whole fragment if a primary lands in an excluded (e.g. rRNA) region,
        // so contaminating reads are never lifted or passed to dedup.
        if(mExcludedRegions != null && mExcludedRegions.fragmentExcluded(group))
        {
            mExcludedReads += group.size();
            return;
        }

        mProcessor.processNameGroup(group, new LiftedMateInfoCache(), this::write);
    }

    public long excludedReads() { return mExcludedReads; }

    private void write(final SAMRecord record, final LiftBackResult result)
    {
        mShardWriter.addAlignment(record);
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
