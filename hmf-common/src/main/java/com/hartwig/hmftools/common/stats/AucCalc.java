package com.hartwig.hmftools.common.stats;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.VectorUtils;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public final class AucCalc
{
    private static final Logger LOGGER = LogManager.getLogger(AucCalc.class);

    public static double calcPercentilesAuc(final List<AucData> dataList, final Level logLevel)
    {
        if(dataList.isEmpty() || dataList.stream().anyMatch(x -> !x.IsPercentile))
            return 0;

        double positiveRankSum = dataList.stream().filter(x -> x.IsPositive).mapToDouble(x -> x.Value).sum();
        long posCount = dataList.stream().filter(x -> x.IsPositive).count();
        return 1 - positiveRankSum / posCount;
    }

    public static double calcScoresAuc(final List<AucData> dataList)
    {
        return calcScoresAuc(dataList, Level.DEBUG);
    }

    public static double calcScoresAuc(final List<AucData> dataList, final Level logLevel)
    {
        if(dataList.isEmpty() || dataList.stream().anyMatch(x -> x.IsPercentile))
            return 0;

        List<Double> posValues = Lists.newArrayList();
        List<Double> negValues = Lists.newArrayList();

        for(AucData data : dataList)
        {
            if(data.IsPositive)
                VectorUtils.optimisedAdd(posValues, data.Value, true);
            else
                VectorUtils.optimisedAdd(negValues, data.Value, true);
        }

        long concordantPairs = 0;
        long discordantPairs = 0;
        long tiedPairs = 0;

        int negIndex = 0;
        int negCount = negValues.size();
        int posCount = posValues.size();

        for(int posIndex = 0; posIndex < posValues.size(); ++posIndex)
        {
            double posValue = posValues.get(posIndex);
            int remainingPos = posCount - posIndex;

            while(negIndex < negValues.size())
            {
                double negValue = negValues.get(negIndex);

                if(Doubles.equal(posValue, negValue))
                    ++tiedPairs;
                else if(negValue < posValue)
                    discordantPairs += remainingPos;
                else
                    break;

                ++negIndex;
            }

            int remainingNeg = negCount - negIndex;
            concordantPairs += remainingNeg;
        }

        long totalPairs = (long)posValues.size() * (long)negValues.size();
        double totalPairsInv = 1.0 / totalPairs;
        double percConcordant = concordantPairs * totalPairsInv;
        double percDiscordant = discordantPairs * totalPairsInv;
        double percTied = tiedPairs * totalPairsInv;

        // Area under curve (AUC) = (Percent Concordant + 0.5 * Percent Tied)/100
        double auc = (percConcordant + 0.5 * percTied);

        if(auc < 0.5)
            auc = 1 - auc;

        LOGGER.log(logLevel, String.format("inputs(pos=%d neg=%d) perc(con=%.4f dis=%.4f tied=%.4f) auc(%.4f)",
                posValues.size(), negValues.size(), percConcordant, percDiscordant, percTied, auc));

        // runSlowMethod(posValues, negValues);

        return auc;
    }

    private static void runSlowMethod(final List<Double> posValues, final List<Double> negValues)
    {
        int concordantPairs = 0;
        int discordantPairs = 0;
        int tiedPairs = 0;

        for(Double posValue : posValues)
        {
            for(Double negValue : negValues)
            {
                if(Doubles.equal(posValue, negValue))
                    ++tiedPairs;
                else if(posValue > negValue)
                    ++concordantPairs;
                else
                    ++discordantPairs;
            }
        }

        double totalPairs = posValues.size() * negValues.size();
        double percConcordant = 100 * concordantPairs / totalPairs;
        double percDiscordant = 100 * discordantPairs / totalPairs;
        double percTied = 100 * tiedPairs / totalPairs;

        // Area under curve (AUC) = (Percent Concordant + 0.5 * Percent Tied)/100
        double auc = (percConcordant + 0.5 * percTied) / 100;

        if(auc < 0.5)
            auc = 1 - auc;

        LOGGER.debug(String.format("slow perc(con=%.4f dis=%.4f tied=%.4f) auc(%.4f)", percConcordant, percDiscordant, percTied, auc));
    }
}
