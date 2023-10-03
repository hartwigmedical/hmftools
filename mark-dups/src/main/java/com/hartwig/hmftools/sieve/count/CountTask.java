package com.hartwig.hmftools.sieve.count;

import static com.hartwig.hmftools.sieve.annotate.AnnotateConfig.MD_LOGGER;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.PerformanceCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CountTask implements Callable
{
    private final CountConfig mConfig;
    private final boolean mDoBucketCounts;
    private final SAMSequenceRecord mSeq;
    private final Vector<CountResult> mResults;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private final PerformanceCounter mPerfCounter;
    private long mReadCounter;

    final private List<Map<Integer, Long>> mPrimaryBucketCounts;
    final private List<Map<Integer, Long>> mSuppBucketCounts;
    private long mPrimaryReadCounter;
    private long mSuppReadCounter;

    public CountTask(final CountConfig config, boolean doBucketCounts, final SAMSequenceRecord seq, final Vector<CountResult> results)
    {
        mConfig = config;
        mDoBucketCounts = doBucketCounts;
        mSeq = seq;
        mResults = results;

        mBamSlicer = new BamSlicer(0, false, true, false);
        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenome)).open(new File(mConfig.BamFile));

        mPerfCounter = new PerformanceCounter(mSeq.getContig());
        mReadCounter = 0;

        mPrimaryReadCounter = 0;
        mSuppReadCounter = 0;

        if(!doBucketCounts)
        {
            mPrimaryBucketCounts = null;
            mSuppBucketCounts = null;
            return;
        }

        final int bucketCount = (mSeq.getSequenceLength() % config.BucketSize == 0)
                ? mSeq.getSequenceLength() / config.BucketSize
                : mSeq.getSequenceLength() / config.BucketSize + 1;
        mPrimaryBucketCounts = new ArrayList<>();
        while(mPrimaryBucketCounts.size() < bucketCount)
        {
            mPrimaryBucketCounts.add(new HashMap<>());
        }

        mSuppBucketCounts = new ArrayList<>();
        while(mSuppBucketCounts.size() < bucketCount)
        {
            mSuppBucketCounts.add(new HashMap<>());
        }
    }

    @Override
    public Long call()
    {
        MD_LOGGER.info("CountTask started for contig {}", mSeq.getContig());

        final ChrBaseRegion chrBaseRegion = new ChrBaseRegion(mSeq.getContig(), 1, mSeq.getSequenceLength());
        mPerfCounter.start();
        mBamSlicer.slice(mSamReader, chrBaseRegion, this::processSamRecord);
        mPerfCounter.stop();
        mResults.add(new CountResult(mSeq, mPrimaryReadCounter, mPrimaryBucketCounts, mSuppReadCounter, mSuppBucketCounts));

        mPerfCounter.logStats();
        MD_LOGGER.info("CountTask for contig {} finished, {} reads processed", mSeq.getContig(), mReadCounter);

        return (long) 0;
    }

    private void processSamRecord(final SAMRecord read)
    {
        ++mReadCounter;

        if(read.getReadUnmappedFlag())
        {
            return;
        }

        if(read.getSupplementaryAlignmentFlag())
        {
            mSuppReadCounter++;
        }
        else
        {
            mPrimaryReadCounter++;
        }

        if(!mDoBucketCounts)
        {
            return;
        }

        final int mapQ = read.getMappingQuality();
        final int bucketIndex = (read.getAlignmentStart() - 1) / mConfig.BucketSize;
        final Map<Integer, Long> bucket =
                (read.getSupplementaryAlignmentFlag()) ? mSuppBucketCounts.get(bucketIndex) : mPrimaryBucketCounts.get(bucketIndex);
        if(!bucket.containsKey(mapQ))
        {
            bucket.put(mapQ, (long) 1);
        }
        else
        {
            bucket.put(mapQ, bucket.get(mapQ) + 1);
        }
    }
}
