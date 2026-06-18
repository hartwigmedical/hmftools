package com.hartwig.hmftools.bamtools.checker;

public class CheckParams
{
    public boolean ConvertHardClips;
    public int MinAlignmentScore;

    public static final int DEFAULT_MIN_ALIGNMENT_SCORE = 30;

    public CheckParams()
    {
        ConvertHardClips = false;
        MinAlignmentScore = DEFAULT_MIN_ALIGNMENT_SCORE;
    }
}
