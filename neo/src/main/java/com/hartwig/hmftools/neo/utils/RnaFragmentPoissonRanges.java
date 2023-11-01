package com.hartwig.hmftools.neo.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.neo.score.TpmCalculator.HIGH_PROBABILITY;
import static com.hartwig.hmftools.neo.score.TpmCalculator.LOW_PROBABILITY;

import com.hartwig.hmftools.neo.score.TpmCalculator;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class RnaFragmentPoissonRanges
{
    public static void main(@NotNull final String[] args) throws ParseException
    {
        Configurator.setRootLevel(Level.INFO);

        System.out.println("Fragments,LowMean,HighMean");

        for(int i = 0; i <= 10; ++i)
        {
            double lowValue = TpmCalculator.calcPoissonMean(i, LOW_PROBABILITY);
            double highValue = TpmCalculator.calcPoissonMean(i, HIGH_PROBABILITY);

            System.out.println(format("%d,%.2f,%.2f", i, lowValue, highValue));
            // NE_LOGGER.info(format("fragments(%d) low(%.2f) high(%.2f)", i, lowValue, highValue));
        }
    }
}
