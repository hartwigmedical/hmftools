package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.neo.bind.RandomDistributionTask.generateDistribution;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class BindUtilsTest
{
    @Test
    public void testScoreDistributions()
    {
        String allele = "B4001";
        int peptideLength = 9;

        List<Double> peptideScores = Lists.newArrayList();

        int totalScoreCount = 1000;
        double scoreStart = 10000;

        for(int i = 0; i < totalScoreCount; ++i)
        {
            peptideScores.add(scoreStart - i);
        }

        List<double[]> discreteScoreData = Lists.newArrayList();
        discreteScoreData.add(new double[] {0.001, 0.1});
        discreteScoreData.add(new double[] {0.01, 0.25});
        discreteScoreData.add(new double[] {0.05, 1.0});

        List<ScoreDistributionData> distributionData = generateDistribution(allele, peptideLength, peptideScores, discreteScoreData);

        assertEquals(131, distributionData.size());

        ScoreDistributionData data = distributionData.get(0);
        assertEquals(peptideScores.get(0), data.Score);
        assertEquals(1, data.CumulativeCount);

        data = distributionData.get(distributionData.size() - 1);
        assertEquals(peptideScores.get(peptideScores.size() - 1), data.Score);
        assertEquals(peptideScores.size(), data.CumulativeCount);
    }
}
