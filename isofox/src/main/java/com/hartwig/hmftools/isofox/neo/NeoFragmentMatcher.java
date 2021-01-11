package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;

import com.hartwig.hmftools.isofox.common.ReadRecord;

public class NeoFragmentMatcher
{
    public static final int MIN_BASE_OVERLAP = 10;

    public static NeoFragmentSupport getNeoEpitopeSupport(final NeoEpitopeData neData, int stream, final ReadRecord read)
    {
        NeoFragmentSupport support = new NeoFragmentSupport();

        // reads always go up in position (+ve to -ve orientation) with any soft-clipped bases potentially on another strand
        final int[] codingBaseRange = new int[SE_PAIR];
        final String neoCodingBases = neData.getFullCodingBases(stream, codingBaseRange);

        // reads with soft-clipped bases to another gene will match the coding bases starting midway through the coding bases
        if(neData.singleGene())
        {
            // the read may cross any part of of the coding base region, which will cover all neo-epitope bases (up, down & novel)
            int maxStartPos = max(read.PosStart, codingBaseRange[SE_START]);
            int minEndPos = min(read.PosEnd, codingBaseRange[SE_END]);

            if(minEndPos - maxStartPos < MIN_BASE_OVERLAP - 1)
                return support;

            final String readBases = extractReadBases(read, maxStartPos, minEndPos);





        }

        // expect the read to either fully fall within one of the up or down stream ranges, or be soft-clippedf




        return support;
    }

    private static String extractReadBases(final ReadRecord read, int posStart, int posEnd)
    {
        int readBaseIndex = 0;
        String bases = "";

        if(read.isSoftClipped(SE_START))
            readBaseIndex += read.Cigar.getFirstCigarElement().getLength();

        for(int[] mappedCoords : read.getMappedRegionCoords(false))
        {
            if(posStart > mappedCoords[SE_END])
                continue;

            for(int basePos = mappedCoords[SE_START]; basePos <= mappedCoords[SE_END]; ++basePos)
            {
                if(basePos > posEnd)
                    return bases;

                if(positionWithin(basePos, posStart, posEnd))
                    bases += read.ReadBases.substring(readBaseIndex, readBaseIndex + 1);

                ++readBaseIndex;
            }
        }

        return bases;
    }


}
