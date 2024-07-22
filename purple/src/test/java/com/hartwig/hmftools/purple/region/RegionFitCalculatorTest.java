package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.purple.PurpleConstants.BAF_PNT_5;
import static com.hartwig.hmftools.purple.fitting.RegionFitCalculator.estimateMinMaxBaf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltTestUtils;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.FittingConfig;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.purple.fitting.RegionFitCalculator;

import org.junit.Test;

public class RegionFitCalculatorTest
{
    private final CobaltChromosomes male = CobaltTestUtils.male();
    private final CobaltChromosomes female = CobaltTestUtils.female();

    @Test
    public void testFitYChromosome()
    {
        final GenomeRegion region = GenomeRegions.create("Y", 1, 100);
        assertTrue(RegionFitCalculator.isAllowedRegion(male, region));
        assertFalse(RegionFitCalculator.isAllowedRegion(female, region));
    }

    @Test
    public void testFitXChromosome()
    {
        final GenomeRegion region = GenomeRegions.create("X", 1, 100);
        assertTrue(RegionFitCalculator.isAllowedRegion(male, region));
        assertTrue(RegionFitCalculator.isAllowedRegion(female, region));
    }

    @Test
    public void testBafDeviationCalcs()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        FittingConfig.addConfig(configBuilder);
        FittingConfig fittingConfig = new FittingConfig(configBuilder, false);

        String chromosome = "1";
        MedianRatio medianRatio = new MedianRatio(chromosome, 0.5, 1);
        CobaltChromosomes cobaltChromosomes = new CobaltChromosomes(Lists.newArrayList(medianRatio), true);

        RegionFitCalculator regionFitCalculator = new RegionFitCalculator(cobaltChromosomes, fittingConfig, 100);
        double purity = 0.08;
        double normFactor = 0.986;

        PurityAdjuster purityAdjuster = new PurityAdjuster(purity, normFactor, cobaltChromosomes);

        double baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 2.018, 0.53);
        assertEquals(0.504, baf, 0.001);

        baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 1.5, 0.53);
        assertEquals(0.667, baf, 0.001);

        baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 2.02, 0.51);
        assertEquals(0.504, baf, 0.001);

        baf = regionFitCalculator.bafToMinimiseDeviation(purityAdjuster, chromosome, 1.9, 0.53);
        assertEquals(0.526, baf, 0.001);
    }

    @Test
    public void testMinMaxBafEstimates()
    {
        // inputs are: copy number, min and max BAF bounds and estimate major allele CN
        double maxBaf = 1;
        double minBaf = BAF_PNT_5;

        testMinMaxBafEstimate(new double[] { 1.25, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 1.5, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 1.75, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 2.0, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 2.1, minBaf, maxBaf, 1.1 });
        testMinMaxBafEstimate(new double[] { 2.25, minBaf, maxBaf, 1.25 });
        testMinMaxBafEstimate(new double[] { 2.5, minBaf, maxBaf, 1.5 });
        testMinMaxBafEstimate(new double[] { 2.75, minBaf, maxBaf, 1.75 });
        testMinMaxBafEstimate(new double[] { 3, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 3.25, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 3.5, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 3.75, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 4.0, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 4.25, minBaf, maxBaf, 2.25 });
        testMinMaxBafEstimate(new double[] { 4.5, minBaf, maxBaf, 2.5 });
        testMinMaxBafEstimate(new double[] { 4.75, minBaf, maxBaf, 2.75 });

        // test with more limited ranges for max BAF
        maxBaf = 0.7;

        testMinMaxBafEstimate(new double[] { 1.25, minBaf, maxBaf, -1 });
        testMinMaxBafEstimate(new double[] { 1.5, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 1.75, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 2.0, minBaf, maxBaf, 1.0 });
        testMinMaxBafEstimate(new double[] { 2.1, minBaf, maxBaf, 1.1 });
        testMinMaxBafEstimate(new double[] { 2.25, minBaf, maxBaf, 1.25 });
        testMinMaxBafEstimate(new double[] { 2.5, minBaf, maxBaf, 1.5 });
        testMinMaxBafEstimate(new double[] { 2.75, minBaf, maxBaf, 1.75 });
        testMinMaxBafEstimate(new double[] { 3, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 3.25, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 3.5, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 3.75, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 4.0, minBaf, maxBaf, 2.0 });
        testMinMaxBafEstimate(new double[] { 4.25, minBaf, maxBaf, 2.25 });
        testMinMaxBafEstimate(new double[] { 4.5, minBaf, maxBaf, 2.5 });
        testMinMaxBafEstimate(new double[] { 4.75, minBaf, maxBaf, 2.75 });
    }

    private static void testMinMaxBafEstimate(final double[] testInputs)
    {
        double estimateBaf = testInputs[3] / testInputs[0];
        double calcBaf = estimateMinMaxBaf(testInputs[0], testInputs[1], testInputs[2]);

        if(calcBaf == -1)
            assertEquals(calcBaf, testInputs[3], 0.001);
        else
            assertEquals(estimateBaf, calcBaf, 0.001);
    }
}
