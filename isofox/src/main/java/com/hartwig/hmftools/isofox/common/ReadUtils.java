package com.hartwig.hmftools.isofox.common;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.bam.CigarUtils.getReadIndexFromPosition;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

public final class ReadUtils
{
    public static void trimAdapterBases(final Read read1, final Read read2)
    {
        if(read1.orientation() == read2.orientation())
            return;

        // no overlap
        if(!positionsOverlap(read1.PosStart, read1.PosEnd, read2.PosStart, read2.PosEnd))
            return;

        // determine the soft-clip bases to trim from each end
        // note that differing splicing can mean that alignments are no the true indication of 5' base boundaries
        int readRefPositionStart1 = 0;
        int readRefPositionStart2 = 0;
        for(int[] mappedCoords1 : read1.getMappedRegionCoords())
        {
            for(int[] mappedCoords2 : read2.getMappedRegionCoords())
            {
                if(positionsOverlap(mappedCoords1[0], mappedCoords1[1], mappedCoords2[0], mappedCoords2[1]))
                {
                    readRefPositionStart1 = readRefPositionStart2 = max(mappedCoords1[0], mappedCoords2[0]);
                    break;
                }
            }

            if(readRefPositionStart1 > 0)
                break;
        }

        if(readRefPositionStart1 == 0)
            return;

        int readIndexStart1 = getReadIndexFromPosition(read1.PosStart, read1.cigarElements(), readRefPositionStart1);
        int readIndexStart2 = getReadIndexFromPosition(read2.PosStart, read2.cigarElements(), readRefPositionStart2);
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
}
