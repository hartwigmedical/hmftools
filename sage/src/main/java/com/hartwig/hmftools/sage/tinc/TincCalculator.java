package com.hartwig.hmftools.sage.tinc;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_GERMLINE_MAX_AD;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_MIN_VARIANTS;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class TincCalculator
{
    private static final int TEST_LEVEL_EARLY_EXIT = 5;

    public static double calculate(final List<VariantData> variants, final List<Double> testLevels)
    {
        if(variants.size() < TINC_MIN_VARIANTS)
            return 0;

        double lowestScore = -1;
        double lowestTestLevel = -1;
        int levelsSinceLowest = 0;

        for(double testLevel : testLevels)
        {
            double score = computeContaminationScore(variants, testLevel);

            SG_LOGGER.trace(format("tinc test level(%.3f) score(%.0f)", testLevel, score));

            if(lowestScore < 0 || score < lowestScore)
            {
                lowestScore = score;
                lowestTestLevel = testLevel;
                levelsSinceLowest = 0;
            }
            else if(levelsSinceLowest >= TEST_LEVEL_EARLY_EXIT)
            {
                // early exit from computing higher levels where the score will continue to worsen
                break;
            }
        }

        SG_LOGGER.debug(format("tinc level(%.3f) minPenalty(%.0f)", lowestTestLevel, lowestScore));

        return lowestTestLevel;
    }

    private static double computeContaminationScore(final List<VariantData> variants, double testLevel)
    {
        /*
            For each variant used in fit, calculate expectedTincAf = TumAF * testedTinc
            Simulate a random variable tincAD = Binom(n=GL_DP, p=expectedTincAf) for each variant
            Get a count of variants by tincAD. Compare this to the count of variants by GL_AD
            Fit penalty is the sum of absolute differences between variant counts at each AD
        */

        double fitPenalty = 0;

        for(int testAd = 0; testAd <= TINC_GERMLINE_MAX_AD; ++testAd)
        {
            double testAdProbabilityTotal = 0;
            int testAdVariantCount = 0;

            for(VariantData variant : variants)
            {
                double tumorAf = variant.tumorAf();
                double expectedAf = testLevel * tumorAf;

                BinomialDistribution distribution = new BinomialDistribution(variant.ReferenceDepth, expectedAf);
                double prob = distribution.probability(testAd);

                testAdProbabilityTotal += prob;

                if(variant.ReferenceAltFrags == testAd
                || (testAd == TINC_GERMLINE_MAX_AD && variant.ReferenceAltFrags > TINC_GERMLINE_MAX_AD))
                {
                    ++testAdVariantCount;
                }
            }

            fitPenalty += abs(testAdVariantCount - testAdProbabilityTotal);
        }

        return fitPenalty;
    }

    private static final double TEST_LEVEL_INCREMENT_1 = 0.005;
    private static final double TEST_LEVEL_INCREMENT_2 = 0.01;
    private static final double TEST_LEVEL_INCREMENT_3 = 0.05;

    private static final double TEST_LEVEL_THRESHOLD_1 = 0.03;
    private static final double TEST_LEVEL_THRESHOLD_2 = 0.15;
    private static final double TEST_LEVEL_THRESHOLD_3 = 1.0;

    public static List<Double> populateDefaultLevels()
    {
        List<Double> testLevels = Lists.newArrayList();

        double increment = TEST_LEVEL_INCREMENT_1;
        for(double level = 0; level < TEST_LEVEL_THRESHOLD_1; level = level + increment)
        {
            testLevels.add(roundedLevel(level, increment));
        }

        increment = TEST_LEVEL_INCREMENT_2;
        for(double level = TEST_LEVEL_THRESHOLD_1; level < TEST_LEVEL_THRESHOLD_2; level = level + increment)
        {
            testLevels.add(roundedLevel(level, 0.01));
        }

        increment = TEST_LEVEL_INCREMENT_3;
        for(double level = TEST_LEVEL_THRESHOLD_2; level <= TEST_LEVEL_THRESHOLD_3; level = level + increment)
        {
            testLevels.add(level);
        }

        return testLevels;
    }

    private static double roundedLevel(double level, double increment)
    {
        return (int)Math.round(level / increment) * increment;
    }
}
