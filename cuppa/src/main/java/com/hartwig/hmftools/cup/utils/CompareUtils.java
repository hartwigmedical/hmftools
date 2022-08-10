package com.hartwig.hmftools.cup.utils;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.cuppa.CuppaDataFile;
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

    public static String topRefResult(final SampleResult result)
    {
        List<String> rankedCancerTypes = CuppaDataFile.getRankedCancerTypes(result.CancerTypeValues);
        return rankedCancerTypes.get(0);
    }

    public static List<String> getRankedCancerTypes(final SampleResult result)
    {
        return CuppaDataFile.getRankedCancerTypes(result.CancerTypeValues);
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
