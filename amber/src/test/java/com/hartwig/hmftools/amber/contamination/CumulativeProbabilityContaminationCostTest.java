package com.hartwig.hmftools.amber.contamination;

import java.util.List;

import com.hartwig.hmftools.amber.VafReading;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.junit.Assert;
import org.junit.Test;

public class CumulativeProbabilityContaminationCostTest
{
    @Test
    public void costForASinglePointTest()
    {
        checkCost(0.08316, 100, 10, 0.1);
        checkCost(0.02656, 1000, 100, 0.1);
        checkCost(0.02656, 1000, 100, 0.2);

        checkCost(0.08316, 100, 90, 0.1);
        checkCost(0.02656, 1000, 900, 0.1);

    }

    @Test
    public void singleValuesFromSpreadsheetTest()
    {
        checkCost(0.5, 6658, 1314, 0.1);
        checkCost(0.43839, 6658, 1314, 0.19);
        checkCost(0.47099, 5858, 169, 0.05);
    }

    @Test
    public void costForMultiplePointsTest()
    {
        VafReading reading1 = rdd(100, 20);
        VafReading reading2 = rdd(80, 16);
        VafReading reading3 = rdd(1000, 200);
        double p = 0.2;
        CumulativeProbabilityContaminationCost model = new CumulativeProbabilityContaminationCost(List.of(reading1, reading2, reading3));
        double e1 = new BinomialDistribution(100, p).cumulativeProbability(20);
        double e2 = new BinomialDistribution(80, p).cumulativeProbability(16);
        double e3 = new BinomialDistribution(1000, p).cumulativeProbability(200);
        Assert.assertEquals(e1 + e2 + e3 - 1.5, model.calculate(p).score(), 0.0001);
    }

    private void checkCost(double expected, int readDepth, int altSupport, double contamination)
    {
        VafReading singleReading = rdd(readDepth, altSupport);
        CumulativeProbabilityContaminationCost model = new CumulativeProbabilityContaminationCost(List.of(singleReading));
        Assert.assertEquals(expected, model.calculate(contamination).score(), 0.0001);
    }

    private VafReading rdd(int depth, int altSupport)
    {
        return new VafReading("1", 100, depth, depth - altSupport, altSupport);
    }
}
