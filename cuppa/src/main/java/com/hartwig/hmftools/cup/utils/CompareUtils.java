package com.hartwig.hmftools.cup.utils;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.cup.common.SampleResult;

public final class CompareUtils
{
    public static boolean resultsMatch(final SampleResult first, final SampleResult second)
    {
        return first.Category == second.Category && first.DataType.equals(second.DataType) && first.Result == second.Result;
    }

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
