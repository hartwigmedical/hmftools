package com.hartwig.hmftools.cup.common;

import java.util.List;

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

    public static void recordCssSimilarity(
            final List<SampleSimilarity> topMatches, final String sample, final String otherSampleId, double css, final String type,
            int maxMatches, double matchCssCutoff)
    {
        if(css < matchCssCutoff)
            return;

        if(topMatches.size() < maxMatches)
        {
            topMatches.add(new SampleSimilarity(sample, otherSampleId, type, css));
            return;
        }

        if(css < topMatches.get(topMatches.size() - 1).Score)
            return;

        for(int i = 0; i < topMatches.size(); ++i)
        {
            if(css > topMatches.get(i).Score)
            {
                topMatches.add(i, new SampleSimilarity(sample, otherSampleId, type, css));
                topMatches.remove(topMatches.size() - 1);
                return;
            }
        }
    }
}
