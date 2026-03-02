package com.hartwig.hmftools.amber.contamination;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.amber.VafReading;

import org.junit.Assert;
import org.junit.Test;

public class LogLikelihoodContaminationCostTest
{
    @Test
    public void altFrequency0FixedDepthTest()
    {
        double c0 = cost(vafReading(100, 0, 0.0), 0.001);
        double c1 = cost(vafReading(100, 1, 0.0), 0.001);
        double c2 = cost(vafReading(100, 2, 0.0), 0.001);
        assertTrue(c0 > c1);
        assertTrue(c1 > c2);
    }

    @Test
    public void altFrequency0FixedAltSupportTest()
    {
        double c0 = cost(vafReading(100, 1, 0.0), 0.001);
        double c1 = cost(vafReading(10, 1, 0.0), 0.001);
        double c2 = cost(vafReading(1, 1, 0.0), 0.001);
        assertTrue(c0 > c1);
        assertTrue(c2 > c1);
    }

    @Test
    public void altFrequency0VariableContaminationLevelTest()
    {
        double c0 = cost(vafReading(100, 1, 0.0), 0.1);
        double c1 = cost(vafReading(100, 1, 0.0), 0.01);
        double c2 = cost(vafReading(100, 1, 0.0), 0.001);
        Assert.assertEquals(c0, c1, 0.00001);
        Assert.assertEquals(c1, c2, 0.00001);
    }

    @Test
    public void altFrequencyLowVariableContaminationLevelTest()
    {
        double c0 = cost(vafReading(100, 1, 0.1), 0.001);
        double c1 = cost(vafReading(100, 1, 0.1), 0.01);
        double c2 = cost(vafReading(100, 1, 0.1), 0.1);
        assertTrue(c0 < c1);
        assertTrue(c2 < c1);
    }

    @Test
    public void altFrequency50PercentTest()
    {
        double c0 = cost(vafReading(100, 1, 0.5), 0.001);
        double c1 = cost(vafReading(100, 1, 0.5), 0.01);
        double c2 = cost(vafReading(100, 1, 0.5), 0.1);
        assertTrue(c0 < c1);
        assertTrue(c2 < c1);
    }

    @Test
    public void costMaximalForExpectedDepth()
    {
        List<Double> costs = new ArrayList<>();
        for(int a = 1; a < 20; a++)
        {
            costs.add(cost(vafReading(100, a, 0.5), 0.1));
        }
        Assert.assertEquals(costs.get(4), costs.stream().max(Double::compareTo).get());
    }

    @Test
    public void scoreDecreasesForFixedRatioAsDepthIncreasesTest()
    {
        double c0 = cost(vafReading(100, 4, 0.5), 0.01);
        double c1 = cost(vafReading(1000, 40, 0.5), 0.01);
        double c2 = cost(vafReading(10000, 400, 0.5), 0.01);
        assertTrue(c0 > c1);
        assertTrue(c1 > c2);
    }

    @Test
    public void handleDataWhereSampleIsHomozygousAlt()
    {
        assertEquals(cost(vafReading(100, 96, 0.5), 0.01), cost(vafReading(100, 4, 0.5), 0.01), 0.0001);
        assertEquals(cost(vafReading(100, 99, 0.6), 0.01), cost(vafReading(100, 1, 0.6), 0.01), 0.0001);
        assertEquals(cost(vafReading(100, 90, 0.4), 0.1), cost(vafReading(100, 10, 0.4), 0.1), 0.0001);
    }

    private VafReading vafReading(int depth, int altSupport, double altFrequency)
    {
        return new VafReading("1", 100, depth, depth - altSupport, altSupport, altFrequency);
    }

    private double cost(VafReading point, double contaminationLevel)
    {
        return new LogLikelihoodContaminationCost(List.of(point)).calculate(contaminationLevel).score();
    }
}