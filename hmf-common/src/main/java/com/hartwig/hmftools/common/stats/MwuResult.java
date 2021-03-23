package com.hartwig.hmftools.common.stats;

public class MwuResult
{
    public double PValue; // probability that cohort 1 is not different from cohort 2
    public int RankSum1; // sum of ranks for cohort 1
    public int Count1; // count of values in cohort 1
    public boolean Valid;

    public MwuResult(boolean valid)
    {
        PValue = 0;
        RankSum1 = 0;
        Count1 = 0;
        Valid = valid;
    }

    // probability that cohort 1 is not less than cohort 2
    public double pValueLess() { return PValue * 0.5; }

    public double pValueGreater() { return 1 - pValueLess(); }

    public double rankAverage1() { return Valid ? RankSum1 / (double)Count1 : 0; }

    public double rankAverage2(int totalValuesCount)
    {
        if(!Valid)
            return 0;

        double rankTotal = totalValuesCount * (totalValuesCount + 1) * 0.5;
        double rankSum2 = rankTotal - RankSum1;
        return rankSum2 / (double)(totalValuesCount - Count1);
    }
}
