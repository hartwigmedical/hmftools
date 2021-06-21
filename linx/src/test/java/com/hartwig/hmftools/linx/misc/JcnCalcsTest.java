package com.hartwig.hmftools.linx.misc;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.cn.CnJcnCalcs.READ_COUNT_PROB_HIGH;
import static com.hartwig.hmftools.linx.cn.CnJcnCalcs.READ_COUNT_PROB_LOW;
import static com.hartwig.hmftools.linx.cn.CnJcnCalcs.calcPoissonObservedGivenProb;

import static junit.framework.TestCase.assertTrue;

import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class JcnCalcsTest
{
    @Test
    public void testExpectedReadCount()
    {
        final List<int[]> probValues = Lists.newArrayList();

        for(int i = 0; i <= 2000; ++i)
        {
            int probCountLow = calcPoissonObservedGivenProb(i, READ_COUNT_PROB_LOW);
            int probCountHigh = calcPoissonObservedGivenProb(i, READ_COUNT_PROB_HIGH);
            // LNX_LOGGER.info("poison prob: expected({}) low({}) high({})", i, probCountLow, probCountHigh);
            probValues.add(i, new int[] {probCountLow, probCountHigh});
        }

        for(int j = 0; j <= 10; ++j)
        {
            for(int i = 0; i <= 2000; ++i)
            {
                int probCountLow = calcPoissonObservedGivenProb(i, READ_COUNT_PROB_LOW);
                int probCountHigh = calcPoissonObservedGivenProb(i, READ_COUNT_PROB_HIGH);

                final int[] existingCounts = probValues.get(i);

                assertTrue(existingCounts[0] == probCountLow);
                assertTrue(existingCounts[1] == probCountHigh);
            }
        }
    }

}
