package com.hartwig.hmftools.bamtools.compare;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.util.List;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.google.common.collect.Queues;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;

public class PartitionReader
{
    private final CompareConfig mConfig;
    private final ChrBaseRegion mRegion;

    private final SamReader mRefSamReader;
    private final SamReader mNewSamReader;
    private final BamSlicer mBamSlicer;

    private final Queue<List<SAMRecord>> mRefReads;
    private List<SAMRecord> mCurrentRefPosReads;
    private final ReadWriter mReadWriter;
    private boolean mLogReadIds;

    private int mRefReadCount;
    private int mNewReadCount;
    private int mDiffCount;

    public PartitionReader(
            final ChrBaseRegion region, final CompareConfig config, final SamReader refSamReader, final SamReader newSamReader,
            final ReadWriter readWriter)
    {
        mConfig = config;
        mRegion = region;
        mReadWriter = readWriter;

        mRefSamReader = refSamReader;
        mNewSamReader = newSamReader;
        mBamSlicer = new BamSlicer(0, true, true, true);
        mBamSlicer.setKeepUnmapped();

        mRefReads = Queues.newArrayDeque();
        mCurrentRefPosReads = null;

        mRefReadCount = 0;
        mNewReadCount = 0;
        mDiffCount = 0;
        mLogReadIds = !mConfig.LogReadIds.isEmpty();
    }

    public void run()
    {
        BT_LOGGER.debug("processing region({})", mRegion);

        mBamSlicer.slice(mRefSamReader, Lists.newArrayList(mRegion), this::processRefRecord);

        mBamSlicer.slice(mNewSamReader, Lists.newArrayList(mRegion), this::processNewRecord);

        while(!mRefReads.isEmpty())
        {
            // write remaining records
            List<SAMRecord> refReads = mRefReads.remove();
            refReads.forEach(x -> mReadWriter.writeComparison(x, REF_ONLY, ""));
            mDiffCount += refReads.size();
        }

        BT_LOGGER.debug("region({}) complete: refReads({}) newReads({}) diff({})",
                mRegion, mRefReadCount, mNewReadCount, mDiffCount);
    }

    private void processRefRecord(final SAMRecord refRead)
    {
        if(!mRegion.containsPosition(refRead.getAlignmentStart()))
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(refRead.getReadName()))
        {
            BT_LOGGER.debug("specific readId({})", refRead.getReadName());
        }

        if(excludeRead(refRead))
            return;

        ++mRefReadCount;

        if(mCurrentRefPosReads == null || mCurrentRefPosReads.get(0).getAlignmentStart() != refRead.getAlignmentStart())
        {
            mCurrentRefPosReads = Lists.newArrayList(refRead);
            mRefReads.add(mCurrentRefPosReads);
        }
        else
        {
            mCurrentRefPosReads.add(refRead);
        }
    }

    private void processNewRecord(final SAMRecord newRead)
    {
        if(!mRegion.containsPosition(newRead.getAlignmentStart()))
            return;

        if(mLogReadIds && mConfig.LogReadIds.contains(newRead.getReadName()))
        {
            BT_LOGGER.debug("specific readId({})", newRead.getReadName());
        }

        if(excludeRead(newRead))
            return;

        ++mNewReadCount;

        if(mRefReads.isEmpty())
        {
            mReadWriter.writeComparison(newRead, NEW_ONLY, "");
            ++mDiffCount;
            return;
        }

        List<SAMRecord> refReads = mRefReads.peek();
        int refPosStart = refReads.get(0).getAlignmentStart();

        // purge any earlier reads
        if(refPosStart < newRead.getAlignmentStart())
        {
            while(!mRefReads.isEmpty())
            {
                refReads = mRefReads.peek();
                refPosStart = refReads.get(0).getAlignmentStart();

                if(refPosStart >= newRead.getAlignmentStart())
                    break;

                refReads.forEach(x -> mReadWriter.writeComparison(x, REF_ONLY, ""));
                mDiffCount += refReads.size();
                mRefReads.remove();
            }
        }

        // write if before the current ref read
        if(newRead.getAlignmentStart() < refPosStart)
        {
            mReadWriter.writeComparison(newRead, NEW_ONLY, "");
            ++mDiffCount;
            return;
        }

        // the positions now match, so look for an exact read match
        int refIndex = 0;
        while(refIndex < refReads.size())
        {
            SAMRecord refRead = refReads.get(refIndex);

            if(readsMatch(refRead, newRead))
            {
                refReads.remove(refIndex);
                if(refReads.isEmpty())
                    mRefReads.remove();

                return;
            }

            ++refIndex;
        }

        // no match
        mReadWriter.writeComparison(newRead, NEW_ONLY, "");
        ++mDiffCount;
    }

    private static boolean readsMatch(final SAMRecord read1, final SAMRecord read2)
    {
        if(!read1.getReadName().equals(read2.getReadName()))
            return false;

        if(read1.getSupplementaryAlignmentFlag() != read2.getSupplementaryAlignmentFlag())
            return false;

        if(read1.getReadUnmappedFlag() != read2.getReadUnmappedFlag())
            return false;

        if(read1.getReadPairedFlag())
        {
            return read1.getFirstOfPairFlag() == read2.getFirstOfPairFlag();
        }
        else
        {
            return true;
        }
    }

    private static boolean excludeRead(final SAMRecord read)
    {
        return read.isSecondaryAlignment() || read.hasAttribute(CONSENSUS_READ_ATTRIBUTE);
    }
}
