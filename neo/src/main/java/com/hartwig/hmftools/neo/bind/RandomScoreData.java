package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.BindData.DELIM;

import java.util.Map;

public class RandomScoreData
{
    public final String Allele;
    public final int PeptideLength;
    public final double ScoreBucket;
    public final double Score;

    public RandomScoreData(final String allele, final int peptideLength, final double scoreBucket, final double score)
    {
        Allele = allele;
        PeptideLength = peptideLength;
        ScoreBucket = scoreBucket;
        Score = score;
    }

    public static RandomScoreData fromCsv(final String data, final Map<String,Integer> fieldsDataMap)
    {
        final String[] items = data.split(DELIM, -1);

        return new RandomScoreData(
                items[fieldsDataMap.get("Allele")], Integer.parseInt(items[fieldsDataMap.get("PeptideLength")]),
                Double.parseDouble(items[fieldsDataMap.get("PeptideLength")]),
                Double.parseDouble(items[fieldsDataMap.get("Score")]));
    }

}
