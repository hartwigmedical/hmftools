package com.hartwig.hmftools.common.redux;

import com.hartwig.hmftools.common.sequencing.ConsensusType;

public class JitterModelParams
{
    public final String RepeatUnit;
    public final com.hartwig.hmftools.common.sequencing.ConsensusType ConsensusType;
    public final double OptimalScaleRepeat4;
    public final double OptimalScaleRepeat5;
    public final double OptimalScaleRepeat6;
    public final double ScaleFitGradient;
    public final double ScaleFitIntercept;
    public final double MicrosatelliteSkew;

    private final int mRepeatUnitLength;

    public static final int MAX_SPECIFIC_LENGTH_UNIT = 2;
    public static final int MIN_REPEAT_COUNT = 4;
    public static final String REPEAT_UNIT_3_PLUS_LABEL = "3-5bp repeat";

    public JitterModelParams(
            final String repeatUnit, final ConsensusType consensusType,
            final double optimalScaleRepeat4, final double optimalScaleRepeat5, final double optimalScaleRepeat6,
            final double scaleFitGradient, final double scaleFitIntercept, final double microsatelliteSkew)
    {
        RepeatUnit = repeatUnit;
        ConsensusType = consensusType;
        OptimalScaleRepeat4 = optimalScaleRepeat4;
        OptimalScaleRepeat5 = optimalScaleRepeat5;
        OptimalScaleRepeat6 = optimalScaleRepeat6;
        ScaleFitGradient = scaleFitGradient;
        ScaleFitIntercept = scaleFitIntercept;
        MicrosatelliteSkew = microsatelliteSkew;

        if(RepeatUnit.contains(REPEAT_DELIM))
        {
            String[] repeats = RepeatUnit.split(REPEAT_DELIM, -1);
            mRepeatUnitLength = repeats[0].length();
        }
        else
        {
            mRepeatUnitLength = MAX_SPECIFIC_LENGTH_UNIT + 1;
        }
    }

    private static final String REPEAT_DELIM = "/";

    public int repeatUnitLength() { return mRepeatUnitLength; }
    public boolean aboveSpecificLength() { return mRepeatUnitLength > MAX_SPECIFIC_LENGTH_UNIT; }
}
