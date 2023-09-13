package com.hartwig.hmftools.bamtools.compare;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.NEW_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.REF_ONLY;
import static com.hartwig.hmftools.bamtools.compare.MismatchType.VALUE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.SUPPLEMENTARY_ATTRIBUTE;

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
            refReads.forEach(x -> mReadWriter.writeComparison(x, REF_ONLY, null));
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
            mReadWriter.writeComparison(newRead, NEW_ONLY, null);
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

                refReads.forEach(x -> mReadWriter.writeComparison(x, REF_ONLY, null));
                mDiffCount += refReads.size();
                mRefReads.remove();
            }
        }

        // write if before the current ref read
        if(newRead.getAlignmentStart() < refPosStart)
        {
            mReadWriter.writeComparison(newRead, NEW_ONLY, null);
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

                checkReadDetails(refRead, newRead);
                return;
            }

            ++refIndex;
        }

        // no match
        mReadWriter.writeComparison(newRead, NEW_ONLY, null);
        ++mDiffCount;
    }

    private static final List<String> KEY_ATTRIBUTES = List.of(SUPPLEMENTARY_ATTRIBUTE, MATE_CIGAR_ATTRIBUTE);

    private void checkReadDetails(final SAMRecord read1, final SAMRecord read2)
    {
        if(read1.getInferredInsertSize() == read2.getInferredInsertSize()
        && read1.getMappingQuality() == read2.getMappingQuality()
        && read1.getFlags() == read2.getFlags()
        && read1.getCigarString().equals(read2.getCigarString()))
        {
            // assume most reads match to avoid creating a array for the diffs
            return;
        }

        List<String> diffs = Lists.newArrayListWithExpectedSize(4);;

        if(read1.getInferredInsertSize() != read2.getInferredInsertSize())
        {
            diffs.add(format("insertSize(%d/%d)", read1.getInferredInsertSize(), read2.getInferredInsertSize()));
        }

        if(read1.getMappingQuality() != read2.getMappingQuality())
        {
            diffs.add(format("mapQuality(%d/%d)", read1.getMappingQuality(), read2.getMappingQuality()));
        }

        if(!read1.getCigarString().equals(read2.getCigarString()))
        {
            diffs.add(format("cigar(%s/%s)", read1.getCigarString(), read2.getCigarString()));
        }

        // check key attributes:
        for(String attribute : KEY_ATTRIBUTES)
        {
            String readAttr1 = read1.getStringAttribute(attribute);
            String readAttr2 = read2.getStringAttribute(attribute);

            if(readAttr1 == null && readAttr2 == null)
                continue;

            if(readAttr1 != null && readAttr2 != null)
            {
                if(!readAttr1.equals(readAttr2))
                {
                    diffs.add(format("attrib_%s(%s/%s)", readAttr1, readAttr2));
                }
            }
            else if(readAttr1 == null && readAttr2 != null)
            {
                diffs.add(format("attrib_%s(missing/%s)", readAttr2));
            }
            else if(readAttr1 != null && readAttr2 == null)
            {
                diffs.add(format("attrib_%s(%s/missing)", readAttr1));
            }
        }

        if(diffs.isEmpty())
            return;

        mReadWriter.writeComparison(read1, VALUE, diffs);
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
