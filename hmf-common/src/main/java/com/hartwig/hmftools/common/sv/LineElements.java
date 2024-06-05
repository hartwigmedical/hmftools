package com.hartwig.hmftools.common.sv;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import com.hartwig.hmftools.common.genome.region.Orientation;

public final class LineElements
{
    public static final String POLY_A_HOMOLOGY = "AAAAAAA";
    public static final String POLY_T_HOMOLOGY = "TTTTTTT";

    public static final int LINE_POLY_AT_TEST_LEN = 18;
    public static final int LINE_POLY_AT_REQ = 16;

    public static boolean isMobileLineElement(final Orientation orientation, final String insertSequence)
    {
        return isMobileLineElement(orientation.asByte(), insertSequence);
    }

    public static boolean isMobileLineElement(final byte orientation, final String insertSequence)
    {
        int insSeqLength = insertSequence.length();
        if(insSeqLength < LINE_POLY_AT_REQ)
            return false;

        final char polyATChar = orientation == POS_ORIENT ? 'T' : 'A';

        int testLength = min(LINE_POLY_AT_TEST_LEN, insSeqLength);
        int allowedNonRequiredChars = testLength - LINE_POLY_AT_REQ;

        for(int i = 0; i < testLength; ++i)
        {
            if(orientation == POS_ORIENT)
            {
                if(insertSequence.charAt(i) != polyATChar)
                    --allowedNonRequiredChars;
            }
            else
            {
                if(insertSequence.charAt(insSeqLength - i - 1) != polyATChar)
                    --allowedNonRequiredChars;
            }

            if(allowedNonRequiredChars < 0)
                return false;
        }

        return true;
    }
}
