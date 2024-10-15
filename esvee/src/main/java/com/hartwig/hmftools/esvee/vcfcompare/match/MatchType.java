package com.hartwig.hmftools.esvee.vcfcompare.match;

public enum MatchType
{
    EXACT_MATCH,
    COORDS_ONLY,
    APPROX_MATCH,
    NO_MATCH;

    public static MatchFunctions.MatchFunction getMatcher(MatchType matchType)
    {
        switch(matchType)
        {
            case EXACT_MATCH: return new MatchFunctions.ExactMatcher();
            case COORDS_ONLY: return new MatchFunctions.CoordsOnlyMatcher();
            case APPROX_MATCH: return new MatchFunctions.ApproxMatcher();
            default: throw new IllegalArgumentException("Invalid match type: " + matchType);
        }
    }
}
