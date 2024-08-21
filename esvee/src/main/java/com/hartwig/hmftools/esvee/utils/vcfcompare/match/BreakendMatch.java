package com.hartwig.hmftools.esvee.utils.vcfcompare.match;

import com.hartwig.hmftools.esvee.utils.vcfcompare.VariantBreakend;

public class BreakendMatch
{
    public VariantBreakend OldBreakend;
    public VariantBreakend NewBreakend;
    public MatchType Type;

    public BreakendMatch(VariantBreakend oldBreakend, VariantBreakend newBreakend, MatchType type)
    {
        OldBreakend = oldBreakend;
        NewBreakend = newBreakend;
        Type = type;
    }
}
