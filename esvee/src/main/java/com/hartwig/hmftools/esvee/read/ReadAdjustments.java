package com.hartwig.hmftools.esvee.read;

import static com.hartwig.hmftools.esvee.SvConstants.NORMALISE_POLY_G_LENGTH;
import static com.hartwig.hmftools.esvee.common.BaseType.G;
import static com.hartwig.hmftools.esvee.common.BaseType.C;

public class ReadAdjustments
{

    public ReadAdjustments()
    {

    }

    public static boolean trimPolyGSequences(final Read read)
    {
        if(read.isUnmapped() || !read.isPairedRead())
            return false;

        int trailingGCount = 0;

        if(read.positiveStrand())
        {
            for(int i = read.basesLength() - 1; i >= 0; --i)
            {
                if(read.getBases()[i] == G.Byte)
                    trailingGCount++;
                else
                    break;
            }
        }
        else
        {
            for(int i = 0; i < read.basesLength(); ++i)
            {
                if(read.getBases()[i] == C.Byte)
                    trailingGCount++;
                else
                    break;
            }
        }

        if(trailingGCount < NORMALISE_POLY_G_LENGTH)
            return false;

        read.trimBases(trailingGCount, read.negativeStrand());

        return true;
    }


}
