package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.max;

import static com.hartwig.hmftools.bamtools.checker.CheckConfig.LOG_READ_COUNT;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.genome.chromosome.Chromosome.isAltRegionContig;

import static org.apache.logging.log4j.Level.DEBUG;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.logging.log4j.Level;

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

    private final FragmentStats mCurrentStats;
    private final FragmentStats mStats;
    private long mReadCount;
    private long mNextLogReadCount;
    private long mCompleteFragments;
    private final boolean mLogReadIds;

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
        mLogReadIds = !mConfig.LogReadIds.isEmpty();

        mCurrentStats = new FragmentStats();
        mStats = new FragmentStats();
    }

    public FragmentStats stats() { return mStats; }

    public void processPartition(final ChrBaseRegion region)
    {
        Level loglevel = isAltRegionContig(region.Chromosome) ? Level.TRACE : DEBUG;
        BT_LOGGER.log(loglevel, "processing region({})", region);

        mRegion = region;

        mReadCount = 0;
        mCompleteFragments = 0;
        mNextLogReadCount = LOG_READ_COUNT;
        mFragmentMap.clear();

        mCurrentStats.reset();

        mBamSlicer.slice(mSamReader, mRegion, this::processSamRecord);

        handleIncompleteFragments();
        mStats.merge(mCurrentStats);

        if(loglevel == Level.TRACE && mReadCount > 0)
            loglevel = DEBUG;

        BT_LOGGER.log(loglevel, "region({}) complete, reads({}) fragments(complete={} incomplete={})",
                mRegion, mReadCount, mCompleteFragments, mFragmentMap.size());

        mFragmentMap.clear();
    }

    private void processSamRecord(final SAMRecord read)
    {
        if(mLogReadIds && mConfig.LogReadIds.contains(read.getReadName()))
        {
            BT_LOGGER.debug("specific read({})", readToString(read));
        }

        if(!mRegion.containsPosition(read.getAlignmentStart()))
            return;

        ++mReadCount;

        if(mReadCount >= mNextLogReadCount)
        {
            mNextLogReadCount += LOG_READ_COUNT;

            BT_LOGGER.info("region({}) processed reads({}) cached fragments({})", mRegion, mReadCount, mFragmentMap.size());
        }

        if(!read.getReadPairedFlag() || read.isSecondaryAlignment())
        {
            writeRecord(read);
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

            mCurrentStats.MaxFragmentCount = max(mCurrentStats.MaxFragmentCount, mFragmentMap.size());
        }
        else
        {
            fragment.addRead(read);
        }

        List<SAMRecord> completeReads = fragment.extractCompleteReads();

        if(!completeReads.isEmpty())
        {
            completeReads.forEach(x -> writeRecord(x));

            if(fragment.isComplete())
            {
                mFragmentMap.remove(fragment.readId());
                ++mCompleteFragments;

                ++mCurrentStats.TotalFragments;

                if(fragment.receivedSupplementaryCount() > 0)
                    ++mCurrentStats.FragmentsWithSupplementaries;

                if(fragment.requiredMateCigarFix())
                    ++mCurrentStats.MateCigarFixed;
            }
        }
    }

    private void handleIncompleteFragments()
    {
        List<SAMRecord> completeReads = mFragmentCache.handleIncompleteFragments(mFragmentMap.values());

        completeReads.forEach(x -> writeRecord(x));
    }

    private void writeRecord(final SAMRecord read)
    {
        if(mBamWriter != null)
            mBamWriter.addAlignment(read);
    }
}
