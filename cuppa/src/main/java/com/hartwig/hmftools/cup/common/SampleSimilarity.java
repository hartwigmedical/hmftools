package com.hartwig.hmftools.cup.common;

public class SampleSimilarity
{
    public final String SampleId;
    public final String MatchedSampleId;
    public final String MatchType;
    public final double Score;

    public SampleSimilarity(final String sampleId, final String matchedSampleId, final String matchType, final double score)
    {
        SampleId = sampleId;
        MatchedSampleId = matchedSampleId;
        MatchType = matchType;
        Score = score;
    }
}
