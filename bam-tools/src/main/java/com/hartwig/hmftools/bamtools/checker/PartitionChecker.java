package com.hartwig.hmftools.bamtools.checker;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionChecker
{
    private final CheckConfig mConfig;
    private final SamReader mSamReader;
    private final SAMFileWriter mBamWriter;
    private final FragmentCache mFragmentCache;

    private final BamSlicer mBamSlicer;
    private ChrBaseRegion mRegion;

    private final Map<String,Fragment> mFragmentMap;

    private long mReadCount;
    private long mNextLogReadCount;
    private long mCompleteFragments;

    private static final int LOG_READ_COUNT = 10_000_000;

    public PartitionChecker(
            final CheckConfig config, final FragmentCache fragmentCache, final SamReader samReader, final SAMFileWriter bamWriter)
    {
        mConfig = config;
        mSamReader = samReader;
        mBamWriter = bamWriter;
        mFragmentCache = fragmentCache;

        mRegion = null;

        mBamSlicer = new BamSlicer(0, true, true, true);
        mBamSlicer.setKeepHardClippedSecondaries();
        mBamSlicer.setKeepUnmapped();

        mFragmentMap = Maps.newHashMap();

        mReadCount = 0;
        mCompleteFragments = 0;
        mNextLogReadCount = LOG_READ_COUNT;
    }

    public void processPartition(final ChrBaseRegion region)
    {
        BT_LOGGER.debug("processing region({})", region);

        mRegion = region;

        mReadCount = 0;
        mCompleteFragments = 0;
        mNextLogReadCount = LOG_READ_COUNT;
        mFragmentMap.clear();

        mBamSlicer.slice(mSamReader, mRegion, this::processSamRecord);

        handleIncompleteFragments();

        BT_LOGGER.debug("region({}) complete, reads({}) fragments(complete={} incomplete={})",
                mRegion, mReadCount, mCompleteFragments, mFragmentMap.size());

        mFragmentMap.clear();
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(!mRegion.containsPosition(read.getAlignmentStart()))
            return;

        ++mReadCount;

        if(mReadCount >= mNextLogReadCount)
        {
            mNextLogReadCount += LOG_READ_COUNT;

            BT_LOGGER.info("region({}) processed reads({}), cached fragments({})", mRegion, mReadCount, mFragmentMap.size());
        }

        if(!read.getReadPairedFlag())
        {
            mBamWriter.addAlignment(read);
            return;
        }

        // check if the mate is in this partition or not
        checkOrCacheRead(read);
    }

    private void checkOrCacheRead(final SAMRecord read)
    {
        Fragment fragment = mFragmentMap.get(read.getReadName());

        if(fragment == null)
        {
            fragment = new Fragment(read);
            mFragmentMap.put(read.getReadName(), fragment);
        }
        else
        {
            fragment.addRead(read);
        }

        List<SAMRecord> completeReads = fragment.extractCompleteReads();

        if(!completeReads.isEmpty())
        {
            completeReads.forEach(x -> mBamWriter.addAlignment(x));

            if(fragment.isComplete())
            {
                mFragmentMap.remove(fragment.readId());
                ++mCompleteFragments;
            }
        }
    }

    private void handleIncompleteFragments()
    {
        List<SAMRecord> completeReads = mFragmentCache.handleIncompleteFragments(mFragmentMap.values());

        completeReads.forEach(x -> mBamWriter.addAlignment(x));
    }
}
