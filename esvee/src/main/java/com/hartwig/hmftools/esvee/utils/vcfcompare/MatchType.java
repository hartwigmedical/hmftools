package com.hartwig.hmftools.esvee.utils.vcfcompare;

enum MatchType
{
    EXACT_MATCH,
    COORDS_ONLY,
    APPROX_MATCH,
    NO_MATCH;

    public static BreakendMatchers.Matcher getMatcher(MatchType matchType)
    {
        switch(matchType)
        {
            case EXACT_MATCH: return new BreakendMatchers.ExactMatcher();
            case COORDS_ONLY: return new BreakendMatchers.CoordsOnlyMatcher();
            case APPROX_MATCH: return new BreakendMatchers.ApproxMatcher();
            default: throw new IllegalArgumentException("Invalid match type: " + matchType);
        }
    }
}
