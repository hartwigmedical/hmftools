package com.hartwig.hmftools.common.sv;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;

import com.hartwig.hmftools.common.genome.region.Orientation;

public final class LineElements
{
    public static final String POLY_A_HOMOLOGY = "AAAAAAA";
    public static final String POLY_T_HOMOLOGY = "TTTTTTT";

    public static final int LINE_POLY_AT_TEST_LEN = 18;
    public static final int LINE_POLY_AT_REQ = 16;

    public static final char LINE_CHAR_A = 'A';
    public static final char LINE_CHAR_T = 'T';
    public static final byte LINE_BASE_A = (byte)LINE_CHAR_A;
    public static final byte LINE_BASE_T = (byte)LINE_CHAR_T;

    public static boolean isMobileLineElement(final Orientation orientation, final String insertSequence)
    {
        return isMobileLineElement(orientation.asByte(), insertSequence);
    }

    public static boolean isLineBase(final byte base) { return base == LINE_BASE_A || base == LINE_BASE_T; }
    public static boolean isLineChar(final char base) { return base == LINE_CHAR_A || base == LINE_CHAR_T; }

    public static boolean isMobileLineElement(final byte orientation, final String insertSequence)
    {
        int insSeqLength = insertSequence.length();
        if(insSeqLength < LINE_POLY_AT_REQ)
            return false;

        final char polyATChar = orientation == ORIENT_FWD ? LINE_CHAR_T : LINE_CHAR_A;

        int testLength = min(LINE_POLY_AT_TEST_LEN, insSeqLength);
        int allowedNonRequiredChars = testLength - LINE_POLY_AT_REQ;

        for(int i = 0; i < testLength; ++i)
        {
            if(orientation == ORIENT_FWD)
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
