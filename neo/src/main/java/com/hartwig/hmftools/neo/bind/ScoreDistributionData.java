package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE_LEN;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;

import java.util.Map;

public class ScoreDistributionData
{
    public final String Allele;
    public final int PeptideLength;
    public final double ScoreBucket;
    public final double Score;

    public final int BucketCount;
    public final int CumulativeCount;

    public ScoreDistributionData(
            final String allele, final int peptideLength, final double scoreBucket, final double score, int bucketCount, int cumulCount)
    {
        Allele = allele;
        PeptideLength = peptideLength;
        ScoreBucket = scoreBucket;
        Score = score;
        BucketCount = bucketCount;
        CumulativeCount = cumulCount;
    }

    public static ScoreDistributionData fromCsv(final String data, final Map<String,Integer> fieldsDataMap)
    {
        final String[] items = data.split(DELIM, -1);

        return new ScoreDistributionData(
                items[fieldsDataMap.get(FLD_ALLELE)],
                Integer.parseInt(items[fieldsDataMap.get(FLD_PEPTIDE_LEN)]),
                Double.parseDouble(items[fieldsDataMap.get("ScoreBucket")]),
                Double.parseDouble(items[fieldsDataMap.get("Score")]),
                0, 0);
    }

    public String toString() { return String.format("bucket(%f) score(%f)", ScoreBucket, Score); }
}
