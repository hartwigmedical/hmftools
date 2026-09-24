package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_INFO_DELIM;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;

import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.SAMRecord;

public final class ReadUtils
{
    public static void trimAdapterBases(final Read read1, final Read read2)
    {
        if(read1.orientation() == read2.orientation())
            return;

        // no overlap
        if(!positionsOverlap(read1.alignmentStart(), read1.alignmentEnd(), read2.alignmentStart(), read2.alignmentEnd()))
            return;

        // determine the soft-clip bases to trim from each end
        // note that differing splicing can mean that alignments are no the true indication of 5' base boundaries
        int readRefPositionStart1 = 0;
        int readRefPositionStart2 = 0;
        for(BaseRegion mappedCoords1 : read1.getMappedRegionCoords())
        {
            for(BaseRegion mappedCoords2 : read2.getMappedRegionCoords())
            {
                if(positionsOverlap(mappedCoords1.start(), mappedCoords1.end(), mappedCoords2.start(), mappedCoords2.end()))
                {
                    readRefPositionStart1 = readRefPositionStart2 = max(mappedCoords1.start(), mappedCoords2.start());
                    break;
                }
            }

            if(readRefPositionStart1 > 0)
                break;
        }

        if(readRefPositionStart1 == 0)
            return;

        int readIndexStart1 = getReadIndexFromPosition(read1.alignmentStart(), read1.cigarElements(), readRefPositionStart1);
        int readIndexStart2 = getReadIndexFromPosition(read2.alignmentStart(), read2.cigarElements(), readRefPositionStart2);
        int readUpperBaseLength1 = read1.baseLength() - readIndexStart1 - 1;
        int readUpperBaseLength2 = read2.baseLength() - readIndexStart2 - 1;

        int trimLength1 = 0;
        int trimLength2 = 0;

        if(read1.orientation().isForward())
        {
            // trim first read on the upper 3' side and vice versa
            trimLength2 = max(readIndexStart2 - readIndexStart1, 0);
            trimLength1 = max(readUpperBaseLength1 - readUpperBaseLength2, 0);
        }
        else
        {
            trimLength1 = max(readIndexStart1 - readIndexStart2, 0);
            trimLength2 = max(readUpperBaseLength2 - readUpperBaseLength1, 0);
        }

        read1.trimAdapterSoftClipBases(trimLength1);
        read2.trimAdapterSoftClipBases(trimLength2);
    }

    public static ClippedSide clippedSide(final Read read)
    {
        int leftScLength = read.leftClipLength();
        int rightScLength = read.rightClipLength();

        if(leftScLength > 0 && rightScLength > 0)
        {
            return leftScLength >= rightScLength ?
                    new ClippedSide(SE_START, leftScLength, true) : new ClippedSide(SE_END, rightScLength, true);
        }
        else if(leftScLength > 0)
        {
            return new ClippedSide(SE_START, leftScLength, true);
        }
        else
        {
            return new ClippedSide(SE_END, rightScLength, rightScLength > 0);
        }
    }

    public static int consensusDuplicateCount(final SAMRecord record)
    {
        String consensusInfo = record.getStringAttribute(CONSENSUS_READ_ATTRIBUTE);

        if(consensusInfo == null)
            return 0;

        String[] consensusComponents = consensusInfo.split(CONSENSUS_INFO_DELIM, 3);
        int duplicateCount = Integer.parseInt(consensusComponents[0]);

        // return 1 less than the cached count since 2 duplicates means 1 primary and 1 duplicate
        return max(duplicateCount - 1, 1);
    }
}
