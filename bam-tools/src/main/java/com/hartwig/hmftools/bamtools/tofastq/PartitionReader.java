package com.hartwig.hmftools.bamtools.tofastq;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.HashMap;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionReader implements AutoCloseable
{
    private final ToFastqConfig mConfig;

    private final RemoteReadHandler mRemoteReadHandler;

    private final FastqWriterCache mWriterCache;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private ChrBaseRegion mCurrentRegion;
    private final Map<String,SAMRecord> mLocalUnmatchedReads;

    private int mPartitionLocalCount;
    private int mPartitionRemoteCount;
    private final PerformanceCounter mPerfCounter;

    public PartitionReader(final ToFastqConfig config, final FastqWriterCache writerCache, final RemoteReadHandler remoteReadHandler)
    {
        mConfig = config;
        mRemoteReadHandler = remoteReadHandler;

        mWriterCache = writerCache;

        mSamReader = ToFastqUtils.openSamReader(mConfig);

        mBamSlicer = new BamSlicer(0, true, false, false);
        mBamSlicer.setKeepUnmapped();

        mCurrentRegion = null;

        mLocalUnmatchedReads = new HashMap<>();

        mPartitionLocalCount = 0;
        mPartitionRemoteCount = 0;

        mPerfCounter = new PerformanceCounter("PartitionReader");
    }

    private static final int LOG_COUNT = 100;

    @Override
    public void close()
    {
        try
        {
            mSamReader.close();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    public PerformanceCounter perfCounter() { return mPerfCounter; }

    public void processRegion(final ChrBaseRegion region)
    {
        mCurrentRegion = region;

        perfCountersStart();

        mBamSlicer.slice(mSamReader, mCurrentRegion, this::processSamRecord);

        postProcessRegion();
    }

    @VisibleForTesting
    public void postProcessRegion()
    {
        // post-slice clean-up
        processPendingIncompletes();

        perfCountersStop();

        BT_LOGGER.debug("partition({}) complete, reads(local={} remote={})",
                mCurrentRegion, mPartitionLocalCount, mPartitionRemoteCount);

        mPartitionLocalCount = 0;
        mPartitionRemoteCount = 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!mCurrentRegion.containsPosition(read.getAlignmentStart())) // to avoid processing reads from the prior region again
            return;

        if(read.hasAttribute(CONSENSUS_READ_ATTRIBUTE)) // drop any consensus reads
            return;

        // check for hard clip
        if(read.getCigar().containsOperator(CigarOperator.HARD_CLIP))
        {
            BT_LOGGER.error("read: {}, hard clip found, require extra logic to handle", read);
            throw new RuntimeException("hard clip found on read");
        }

        if(!read.getReadPairedFlag())
        {
            ++mPartitionLocalCount;
            mWriterCache.writeUnpairedRead(read);
            return;
        }

        SAMRecord mate = mLocalUnmatchedReads.remove(read.getReadName());

        if(mate != null)
        {
            ++mPartitionLocalCount;
            mWriterCache.writeReadPair(read, mate);
            return;
        }

        boolean hasLocalMate =read.getReferenceIndex().equals(read.getMateReferenceIndex())
                && mCurrentRegion.containsPosition(read.getMateAlignmentStart());

        // cache if local otherwise put into remote pending list
        if(hasLocalMate)
        {
            ++mPartitionLocalCount;
            mLocalUnmatchedReads.put(read.getReadName(), read);
            return;
        }

        // pass any remote read to the remote read handler
        mRemoteReadHandler.cacheRemoteRead(read);
        ++mPartitionRemoteCount;
    }

    private void processPendingIncompletes()
    {
        // any local unmatched read needs to be logged, and passed to remote read handler
        if(!mLocalUnmatchedReads.isEmpty())
        {
            BT_LOGGER.warn("partition({}) has {} unmatched local reads", mCurrentRegion, mLocalUnmatchedReads.size());

            for(SAMRecord read : mLocalUnmatchedReads.values())
            {
                Validate.isTrue(read.getReadPairedFlag());
                BT_LOGGER.error("unmatched local paired read: {}", read);

                // shouldn't happen, pass it to remote read handler
                mRemoteReadHandler.cacheRemoteRead(read);
            }

            mLocalUnmatchedReads.clear();
        }
    }

    private void perfCountersStart()
    {
        if(mConfig.PerfDebug)
            mPerfCounter.start(format("%s", mCurrentRegion));
        else
            mPerfCounter.start();
    }

    private void perfCountersStop()
    {
        mPerfCounter.stop();
    }

    @VisibleForTesting
    public void processRead(final SAMRecord read) { processSamRecord(read); }

    @VisibleForTesting
    public void flushPendingIncompletes() { processPendingIncompletes(); }
}
