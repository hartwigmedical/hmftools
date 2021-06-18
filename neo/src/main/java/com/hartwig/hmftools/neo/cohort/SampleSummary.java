package com.hartwig.hmftools.neo.cohort;

import static com.hartwig.hmftools.neo.cohort.PredictionData.DELIM;
import static com.hartwig.hmftools.neo.cohort.StatusResults.STATUS_MAX;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class SampleSummary
{
    public final int PeptideCount;

    public final StatusResults[] Results;

    public SampleSummary(int peptideCount)
    {
        PeptideCount = peptideCount;
        Results = new StatusResults[STATUS_MAX];

        for(int i = 0; i < STATUS_MAX; ++i)
        {
            Results[i] = new StatusResults();
        }
    }

    public static String header()
    {
        List<String> columns = Lists.newArrayList("AffTotal", "AffLowCount", "AffMedCount", "PresTotal", "PresCount");
        List<String> statusGroups = Lists.newArrayList("Norm", "Tumor", "SimTumor");

        StringJoiner sj = new StringJoiner(DELIM);
        for(String group : statusGroups)
        {
            for(String column : columns)
            {
                sj.add(group + column);
            }
        }

        return sj.toString();
    }

}
