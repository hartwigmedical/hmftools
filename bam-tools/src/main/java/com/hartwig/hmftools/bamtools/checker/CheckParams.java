package com.hartwig.hmftools.bamtools.checker;

public class CheckParams
{
    public boolean ConvertHardClips;
    public int MinSuppAlignmentScore;

    public static final int DEFAULT_MIN_SUPP_ALIGNMENT_SCORE = 30;

    public CheckParams()
    {
        ConvertHardClips = false;
        MinSuppAlignmentScore = DEFAULT_MIN_SUPP_ALIGNMENT_SCORE;
    }
}
