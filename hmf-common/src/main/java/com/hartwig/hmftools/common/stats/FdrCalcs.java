package com.hartwig.hmftools.common.stats;

import java.util.Collections;
import java.util.List;

public class FdrCalcs
{
    public static void calculateFDRs(final List<PValueResult> pValues)
    {
        // rank P-values then calculate a FDR for each of them
        Collections.sort(pValues);
        int testCount = pValues.size();

        for(int i = 0; i < pValues.size(); ++i)
        {
            PValueResult pValue = pValues.get(i);
            pValue.Rank = i + 1;
            pValue.QValue = pValue.PValue * testCount / pValue.Rank;
        }
    }
}
