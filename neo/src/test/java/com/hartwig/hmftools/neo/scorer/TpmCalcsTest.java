package com.hartwig.hmftools.neo.scorer;

import static java.lang.String.format;

import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.scorer.TpmCalculator.HIGH_PROBABILITY;
import static com.hartwig.hmftools.neo.scorer.TpmCalculator.LOW_PROBABILITY;

import static junit.framework.TestCase.assertEquals;

import org.junit.Test;

public class TpmCalcsTest
{
    @Test
    public void testEffectiveTpmCalcs()
    {
        double reqProbLow = 0.05;
        double reqProbHigh = 0.95;

        double lowValue = TpmCalculator.calcPoissonMean(10, reqProbLow);
        double highValue = TpmCalculator.calcPoissonMean(10, reqProbHigh);

        assertEquals(17.0, lowValue, 0.1);
        assertEquals(6.3, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(0, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(0, reqProbHigh);

        assertEquals(3, lowValue, 0.1);
        assertEquals(0, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(1, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(1, reqProbHigh);

        assertEquals(4.7, lowValue, 0.1);
        assertEquals(0.3, highValue, 0.1);

        lowValue = TpmCalculator.calcPoissonMean(2, reqProbLow);
        highValue = TpmCalculator.calcPoissonMean(2, reqProbHigh);

        assertEquals(6.3, lowValue, 0.1);
        assertEquals(0.8, highValue, 0.1);
    }

    @Test
    public void generateCachedValues()
    {
        // Configurator.setRootLevel(Level.INFO);

        for(int i = 0; i <= 10; ++i)
        {
            double lowValue = TpmCalculator.calcPoissonMean(i, LOW_PROBABILITY);
            double highValue = TpmCalculator.calcPoissonMean(i, HIGH_PROBABILITY);

            NE_LOGGER.info(format("fragments(%d) low(%.2f) high(%.2f)", i, lowValue, highValue));
        }
    }

}
