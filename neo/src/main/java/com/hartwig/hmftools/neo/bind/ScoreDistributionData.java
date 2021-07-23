package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindData.DELIM;

import java.util.Map;

public class ScoreDistributionData
{
    public final String Allele;
    public final int PeptideLength;
    public final double ScoreBucket;
    public final double Score;

    public ScoreDistributionData(final String allele, final int peptideLength, final double scoreBucket, final double score)
    {
        Allele = allele;
        PeptideLength = peptideLength;
        ScoreBucket = scoreBucket;
        Score = score;
    }

    public static ScoreDistributionData fromCsv(final String data, final Map<String,Integer> fieldsDataMap)
    {
        final String[] items = data.split(DELIM, -1);

        return new ScoreDistributionData(
                items[fieldsDataMap.get("Allele")], Integer.parseInt(items[fieldsDataMap.get("PeptideLength")]),
                Double.parseDouble(items[fieldsDataMap.get("ScoreBucket")]), Double.parseDouble(items[fieldsDataMap.get("Score")]));
    }

}
