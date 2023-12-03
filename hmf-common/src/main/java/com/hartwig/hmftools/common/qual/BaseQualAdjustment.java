package com.hartwig.hmftools.common.qual;

public class BaseQualAdjustment
{
    public static final int[] STANDARD_BASE_QUALS = { 0, 11, 25, 37 };

    public static final double BASE_QUAL_PERMITTED_DIFF_MIN = 0.1;
    public static final double BASE_QUAL_PERMITTED_DIFF_MAX = 2;

    public static final int BASE_QUAL_LOWER_STEP = 1; // step down from the standard value if difference is within bounds

    public static byte adjustBaseQual(final double baseQual) { return adjustBaseQual(STANDARD_BASE_QUALS, baseQual); }

    public static byte adjustBaseQual(final int[] standardBaseQuals, final double baseQual)
    {
        if(baseQual < standardBaseQuals[1] - BASE_QUAL_PERMITTED_DIFF_MAX)
            return (byte)standardBaseQuals[0];

        // for raw base qual lower than the standard value:
        // - by the min difference or less, use the standard value, eg: 24.9 -> 25
        // - by the max difference or less, use the standard value - step value, eg: 23.1 -> 25 - 1 = 24

        for(int i = standardBaseQuals.length - 1; i > 0; --i)
        {
            if(baseQual > standardBaseQuals[i])
                return (byte)standardBaseQuals[i];

            double diff = standardBaseQuals[i] - baseQual;

            if(diff <= BASE_QUAL_PERMITTED_DIFF_MIN)
                return (byte)standardBaseQuals[i];

            if(diff <= BASE_QUAL_PERMITTED_DIFF_MAX)
                return (byte)(standardBaseQuals[i] - BASE_QUAL_LOWER_STEP);
        }

        return (byte)STANDARD_BASE_QUALS[0];
    }
}
