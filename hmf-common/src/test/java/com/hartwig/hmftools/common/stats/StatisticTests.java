package com.hartwig.hmftools.common.stats;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class StatisticTests
{
    @Test
    public void testStats()
    {
        double[] vals1 = {1,3,5,7,9,11};
        double[] vals2 = {2,4,6,8,10,12};

        MannWhitneyUTest mwu = new MannWhitneyUTest();

        MwuResult result = mwu.calculate(vals1, vals2);

        assertTrue(result.Valid);
        assertEquals(36, result.RankSum1);
        assertEquals(6, result.Count1);
        assertEquals(0.631, result.PValue, 0.001);

        // test again with different size cohorts
        vals1 = new double[] {1,2,3,4,5,6,7,8,12};
        vals2 = new double[] {9,10,11};

        result = mwu.calculate(vals1, vals2);

        assertTrue(result.Valid);
        assertEquals(48, result.RankSum1);
        assertEquals(9, result.Count1);
        assertEquals(0.052, result.PValue, 0.001);

    }
}
