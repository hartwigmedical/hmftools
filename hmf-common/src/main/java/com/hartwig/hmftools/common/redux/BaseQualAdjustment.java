package com.hartwig.hmftools.common.redux;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import com.hartwig.hmftools.common.sequencing.IlluminaBamUtils;
import com.hartwig.hmftools.common.sequencing.SbxBamUtils;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;

public class BaseQualAdjustment
{
    // zero is not handled by some downstream tools
    // this low value also denotes an adjustment (eg by Redux) to indicate that the corresponding base is completely uncertain
    public static final byte BASE_QUAL_MINIMUM = 1;

    public static final byte LOW_BASE_QUAL_THRESHOLD = 26; // quals below this value are considered low-qual

    private static final int[] STANDARD_BASE_QUALS = { BASE_QUAL_MINIMUM, 11, 25, 37 };
    private static final double BASE_QUAL_PERMITTED_DIFF_MAX = 1.5;

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

    public static byte maxQual(final byte qual1, final byte qual2) { return (byte)max(qual1, qual2); }
    public static byte minQual(final byte qual1, final byte qual2) { return (byte)min(qual1, qual2); }

    public static double phredQualToProbability(byte quality)
    {
        return pow(10, -quality / 10.0);
    }

    public static byte probabilityToPhredQualInt(double probability)
    {
        return (byte)round(probabilityToPhredQual(probability));
    }

    public static double probabilityToPhredQual(double probability) { return -10 * log10(probability); }

    public static boolean isLowBaseQual(final byte qual) { return qual < LOW_BASE_QUAL_THRESHOLD; }
    public static boolean aboveLowBaseQual(final byte qual) { return qual >= LOW_BASE_QUAL_THRESHOLD; }

    public static boolean isHighBaseQual(final byte qual, final SequencingType sequencingType)
    {
        if(sequencingType == SequencingType.ILLUMINA)
            return IlluminaBamUtils.isHighBaseQual(qual);
        else if(sequencingType == SequencingType.SBX)
            return SbxBamUtils.isHighBaseQual(qual);
        else
            return UltimaBamUtils.isHighBaseQual(qual);
    }

    public static boolean isMediumBaseQual(final byte qual, final SequencingType sequencingType)
    {
        if(sequencingType == SequencingType.SBX)
            return SbxBamUtils.isMediumBaseQual(qual);
        else
            return false;
    }

    public static boolean isUncertainBaseQual(final byte qual) { return qual <= BASE_QUAL_MINIMUM; }
}
