package com.hartwig.hmftools.bamtools.checker;

import java.util.Set;

import com.google.common.collect.Sets;

public class CheckParams
{
    public boolean ConvertHardClips;
    public int MinAlignmentScore;
    public final Set<String> ValidContigs;

    public static final int DEFAULT_MIN_ALIGNMENT_SCORE = 30;

    public CheckParams()
    {
        ConvertHardClips = false;
        MinAlignmentScore = DEFAULT_MIN_ALIGNMENT_SCORE;
        ValidContigs = Sets.newHashSet();
    }
}
