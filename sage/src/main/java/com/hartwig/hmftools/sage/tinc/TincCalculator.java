package com.hartwig.hmftools.sage.tinc;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class TincCalculator
{
    private final TincConfig mConfig;
    private final List<VariantData> mVariants;
    private final List<Double> mTestLevels;

    protected static final int TINC_GERMLINE_ABQ_MIN = 30;
    protected static final int TINC_MAX_FITTING_VARIANTS = 15_000;

    protected static final double TINC_GERMLINE_DEPTH_LOW = 0.5;
    protected static final double TINC_GERMLINE_DEPTH_HIGH = 1 + TINC_GERMLINE_DEPTH_LOW;

    public TincCalculator(final TincConfig config, final List<VariantData> variants)
    {
        mConfig = config;
        mVariants = variants;

        mTestLevels = Lists.newArrayList();
        populateLevels();
    }

    public void run()
    {
        double lowestScore = -1;
        double lowestTestLevel = 0;

        for(double testLevel : mTestLevels)
        {
            double score = computeContaminationScore(testLevel);

            if(lowestScore < 0 || score < lowestScore)
            {
                lowestScore = score;
                lowestTestLevel = testLevel;
            }
        }
    }

    private double computeContaminationScore(double testLevel)
    {
        /*
            For each variant used in fit, calculate expectedTincAf = TumAF * testedTinc
            Simulate a random variable tincAD = Binom(n=GL_DP, p=expectedTincAf) for each variant
            Get a count of variants by tincAD. Compare this to the count of variants by GL_AD
            Fit penalty is the sum of absolute differences between variant counts at each AD
        */

        double totalProbability = 0;

        /*
        for(VariantData variant : mVariants)
        {
            double expectedAf = testLevel * variant.TumorAltFrags;

            BinomialDistribution distribution = new BinomialDistribution(variant.ReferenceDepth, expectedAf);
            double prob = distribution.cumulativeProbability(variant.ReferenceAltFrags);

            totalProbability += prob;
        }
        */

        return totalProbability;
    }

    private void populateLevels()
    {
        double increment = 0.005;
        for(double level = 0; level < 0.03; level = level + increment)
        {
            mTestLevels.add(roundedLevel(level, increment));
        }

        increment = 0.01;
        for(double level = 0.03; level < 0.15; level = level + increment)
        {
            mTestLevels.add(roundedLevel(level, 0.01));
        }

        increment = 0.05;
        for(double level = 0.15; level <= 1; level = level + increment)
        {
            mTestLevels.add(level);
        }
    }

    private static double roundedLevel(double level, double increment)
    {
        return (int)Math.round(level / increment) * increment;
    }

    public void setTestLevels(final List<Double> levels)
    {
        mTestLevels.clear();
        mTestLevels.addAll(levels);
    }
}
