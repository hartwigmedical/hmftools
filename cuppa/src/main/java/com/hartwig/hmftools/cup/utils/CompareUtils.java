package com.hartwig.hmftools.cup.utils;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.common.SampleResult;

public final class CompareUtils
{
    public static boolean resultsMatch(final SampleResult first, final SampleResult second)
    {
        return first.Category == second.Category && dataTypesMatch(first.DataType, second.DataType) && first.Result == second.Result;
    }

    private static boolean dataTypesMatch(final String dataType1, final String dataType2)
    {
        return dataType1.equals(dataType2); // || (isGenPosDataType(dataType1) && isGenPosDataType(dataType2));
    }

    /*
    private static boolean isGenPosDataType(final String dataType)
    {
        return dataType.equals(GENOMIC_POSITION_SIMILARITY.toString()) || dataType.equals(GENOMIC_POSITION_PAIRWISE.toString());
    }
    */

    public static String topRefResult(final SampleResult result)
    {
        double topValue = 0;
        String topCancerType = "";

        for(Map.Entry<String,Double> entry : result.CancerTypeValues.entrySet())
        {
            if(entry.getValue() > topValue)
            {
                topValue = entry.getValue();
                topCancerType = entry.getKey();
            }
        }

        return topCancerType;
    }

    public static List<String> getRankedCancerTypes(final SampleResult result)
    {
        List<String> cancerTypes = Lists.newArrayList();
        List<Double> scores = Lists.newArrayList();

        for(Map.Entry<String,Double> entry : result.CancerTypeValues.entrySet())
        {
            int index = 0;
            double score = entry.getValue();

            while(index < scores.size())
            {
                if(score > scores.get(index))
                    break;

                ++index;
            }

            scores.add(index, score);
            cancerTypes.add(index, entry.getKey());
        }

        return cancerTypes;
    }


    public static final String EMPTY_RESULTS_CSV = DATA_DELIM + DATA_DELIM;

    public static String resultInfoCsv(final SampleResult result, final String refCancerType)
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);

        String topCancerType = topRefResult(result);
        sj.add(topCancerType);
        sj.add(String.valueOf(result.CancerTypeValues.get(topCancerType)));
        sj.add(String.valueOf(result.CancerTypeValues.get(refCancerType)));

        return sj.toString();
    }

}
