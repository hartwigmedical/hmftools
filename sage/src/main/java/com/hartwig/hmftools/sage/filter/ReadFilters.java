package com.hartwig.hmftools.sage.filter;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.samtools.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.samtools.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.mateNegativeStrand;
import static com.hartwig.hmftools.sage.SageConstants.CHIMERIC_FRAGMENT_LENGTH_MAX;

import com.hartwig.hmftools.common.samtools.SamRecordUtils;

import htsjdk.samtools.SAMRecord;

public final class ReadFilters
{
    public static boolean isChimericRead(final SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            if(record.getMateUnmappedFlag())
                return true;

            // inter-chromosomal
            if(!record.getReferenceName().equals(record.getMateReferenceName()))
                return true;

            // inversion
            if(record.getReadNegativeStrandFlag() == mateNegativeStrand(record))
                return true;
        }

        // or a fragment length outside the expected maximum
        if(abs(record.getInferredInsertSize()) > CHIMERIC_FRAGMENT_LENGTH_MAX)
            return true;

        return false;
    }

    public static final int NO_READ_EDGE_DISTANCE = -1;

    public static int calcDistanceFromReadEdge(final int variantPosition, final SAMRecord record, boolean includeSoftClipBases)
    {
        if(!positionWithin(variantPosition, record.getAlignmentStart(), record.getAlignmentEnd()))
            return NO_READ_EDGE_DISTANCE;

        int distFromStart = variantPosition - record.getAlignmentStart();
        int distFromEnd = record.getAlignmentEnd() - variantPosition;

        if(includeSoftClipBases)
        {
            distFromStart += leftSoftClipLength(record);
            distFromEnd += rightSoftClipLength(record);
        }

        return min(distFromStart, distFromEnd);
    }
}
