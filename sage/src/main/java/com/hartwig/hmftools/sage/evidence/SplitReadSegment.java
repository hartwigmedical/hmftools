package com.hartwig.hmftools.sage.evidence;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;
import static htsjdk.samtools.CigarOperator.S;

import java.util.Arrays;

import com.hartwig.hmftools.common.bam.CigarUtils;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class SplitReadSegment
{
    public final int ReadVarIndex;
    public final int SegmentIndexStart; // from the original read
    public final int SegmentIndexEnd;
    public final byte[] ReadBases;
    public final byte[] ReadQuals;

    public SplitReadSegment(final int readVarIndex, final int segIndexStart, final int segIndexEnd, final byte[] readBases, final byte[] readQuals)
    {
        ReadVarIndex = readVarIndex;
        SegmentIndexStart = segIndexStart;
        SegmentIndexEnd = segIndexEnd;
        ReadBases = readBases;
        ReadQuals = readQuals;
    }

    public int length() { return ReadBases.length; }

    public String toString() { return format("index(%s) seg(%d-%d) baseLength(%d)",
            ReadVarIndex, SegmentIndexStart, SegmentIndexEnd, ReadBases.length); }

    public static SplitReadSegment formSegment(final SAMRecord record, final int position, final int readVarIndex)
    {
        int readIndex = 0;
        int refPosition = record.getAlignmentStart() - CigarUtils.leftSoftClipLength(record);

        int segmentIndexStart = -1;
        int segmentIndexEnd = -1;

        for(CigarElement element : record.getCigar().getCigarElements())
        {
            if(element.getOperator().consumesReferenceBases() || element.getOperator() == S)
            {
                int refEndPosition = refPosition + element.getLength() - 1;

                boolean positionWithinSection = positionWithin(position, refPosition, refEndPosition);

                if(positionWithinSection && (element.getOperator() == N || element.getOperator() == D))
                    return null;

                if(element.getOperator() == N)
                {
                    if(refPosition > position)
                        break;

                    // reset base capture after any split
                    segmentIndexStart = -1;
                    segmentIndexEnd = -1;
                }
                else if(element.getOperator() == M || element.getOperator() == S)
                {
                    if(segmentIndexStart < 0)
                        segmentIndexStart = readIndex;
                }

                refPosition += element.getLength();
            }

            if(element.getOperator().consumesReadBases())
            {
                if(segmentIndexStart >= 0)
                {
                    if(segmentIndexEnd == -1)
                        segmentIndexEnd = segmentIndexStart + element.getLength() - 1;
                    else
                        segmentIndexEnd += element.getLength();
                }

                readIndex += element.getLength();
            }
        }

        if(segmentIndexStart < 0)
            return null;

        byte[] segmentBases = Arrays.copyOfRange(record.getReadBases(), segmentIndexStart, segmentIndexEnd + 1);
        byte[] segmentQuals = Arrays.copyOfRange(record.getBaseQualities(), segmentIndexStart, segmentIndexEnd + 1);

        return new SplitReadSegment(
                readVarIndex - segmentIndexStart, segmentIndexStart, segmentIndexEnd, segmentBases, segmentQuals);
    }
}
