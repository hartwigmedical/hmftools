package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;

import com.hartwig.hmftools.isofox.common.ReadRecord;

public class NeoFragmentMatcher
{
    public static NeoFragmentSupport getNeoEpitopeSupport(final NeoEpitopeData neData, int stream, final ReadRecord read)
    {
        // reads always go up in position (+ve to -ve orientation) with any soft-clipped bases potentially on another strand
        final int[] codingBaseRange = new int[SE_PAIR];
        final String neoCodingBases = neData.getFullCodingBases(stream, codingBaseRange);




        int maxStartPos = max(read.PosStart, neData.Source.CodingBasePositions[FS_DOWN]);


        return new NeoFragmentSupport();
    }


}
