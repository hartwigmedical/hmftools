package com.hartwig.hmftools.common.qual;

public class BaseQualAdjustment
{
    public static final byte BASE_QUAL_MINIMUM = 1; // zero is not handled by some downstream tools
    public static final int[] STANDARD_BASE_QUALS = { BASE_QUAL_MINIMUM, 11, 25, 37 };

    public static final double BASE_QUAL_PERMITTED_DIFF_MAX = 2;

    public static byte adjustBaseQual(final double baseQual) { return adjustBaseQual(STANDARD_BASE_QUALS, baseQual); }

    public static byte adjustBaseQual(final int[] standardBaseQuals, final double baseQual)
    {
        if(baseQual < standardBaseQuals[1] - BASE_QUAL_PERMITTED_DIFF_MAX)
            return (byte)standardBaseQuals[0];

        // for raw base qual lower than the standard value:
        // - by the max difference or less, use the standard value, eg: 24.9 -> 25

        for(int i = standardBaseQuals.length - 1; i > 0; --i)
        {
            if(baseQual > standardBaseQuals[i])
                return (byte)standardBaseQuals[i];

            double diff = standardBaseQuals[i] - baseQual;

            if(diff <= BASE_QUAL_PERMITTED_DIFF_MAX)
                return (byte)standardBaseQuals[i];
        }

        return (byte)STANDARD_BASE_QUALS[0];
    }
}
