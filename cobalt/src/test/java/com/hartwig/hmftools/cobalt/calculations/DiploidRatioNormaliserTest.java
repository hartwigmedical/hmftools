package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

public class DiploidRatioNormaliserTest extends CalculationsTestBase
{
    private static final double EPSILON = 1e-10;
    private int minWindowCoverage = 2;
    private int maxWindowDistance = 2;

    @Test
    public void normaliseTest()
    {
        List<BamRatio> ratios = new ArrayList<>();
        ratios.add(br(1, -1.0));
        ratios.add(br(1001, 0.50));
        ratios.add(br(2001, 0.51));
        ratios.add(br(3001, 0.51));
        ratios.add(br(4001, 0.51));
        ratios.add(br(5001, 0.49));
        ratios.add(br(6001, 0.58));
        ratios.add(br(7001, 0.51));
        ratios.add(br(8001, 0.58));
        ratios.add(br(9001, 0.49));
        ratios.add(br(10001, 0.49));
        ratios.add(br(11001, -1.0));
        ratios.add(br(12001, 0.59));

        DiploidRatioNormaliser normaliser = new DiploidRatioNormaliser(maxWindowDistance, minWindowCoverage);
        ratios.forEach(normaliser::recordRatio);

        normaliser.dataCollectionFinished();
        normaliser.setExpectedRatio(1.0);
        assertEquals(11, normaliser.count());
        assertEquals(0.51, normaliser.median(), 0.0001);

        assertEquals(-1.0, normaliser.normalise(ratios.get(0)), EPSILON);
        assertEquals(0.5/0.51, normaliser.normalise(ratios.get(1)), EPSILON);
        assertEquals(1.0, normaliser.normalise(ratios.get(2)), EPSILON);
        assertEquals(1.0, normaliser.normalise(ratios.get(3)), EPSILON);
        assertEquals(1.0, normaliser.normalise(ratios.get(4)), EPSILON);
        assertEquals(0.49/0.51, normaliser.normalise(ratios.get(5)), EPSILON);
        assertEquals(0.58/0.51, normaliser.normalise(ratios.get(6)), EPSILON);
        assertEquals(1.0, normaliser.normalise(ratios.get(7)), EPSILON);
        assertEquals(0.58/0.51, normaliser.normalise(ratios.get(8)), EPSILON);
        assertEquals(0.49/0.51, normaliser.normalise(ratios.get(9)), EPSILON);
        assertEquals(0.49/0.51, normaliser.normalise(ratios.get(10)), EPSILON);
        assertEquals(-1.0, normaliser.normalise(ratios.get(11)), EPSILON);
        assertEquals(0.59/0.50, normaliser.normalise(ratios.get(12)), EPSILON);
    }

    @Test
    public void testCloseToZero()
    {
        final List<Double> input = Arrays.asList(0.0, 0.0, 0.002, 0.0, 0.0);
        maxWindowDistance = 5;
        minWindowCoverage = 5;
        final List<Double> output = doNormalisation(input, 1.0);//1.0, 5, 5, input).get();
        assertEquals(input.size(), output.size());
        assertRatio(input.get(0), output.get(0), 1);
        assertRatio(input.get(1), output.get(1), 1);
        assertRatio(input.get(2), output.get(2), 1);
        assertRatio(input.get(3), output.get(3), 1);
        assertRatio(input.get(4), output.get(4), 1);
    }

    @Test
    public void testMaxWindowDistance()
    {
        final List<Double> input = Arrays.asList(1.0, 1.5, -1.0, 1.1, 1.2);
        maxWindowDistance = 2;
        minWindowCoverage = 1;
        List<Double> output = doNormalisation(input, 1.15);

        assertEquals(input.size(), output.size());
        assertEquals(1.15 /1.1, output.get(0), EPSILON);
        assertEquals(1.5*1.15/1.15, output.get(1), EPSILON);
        assertEquals(-1.0, output.get(2), EPSILON);
        assertEquals(1.1*1.15/1.15, output.get(3), EPSILON);
        assertEquals(1.2*1.15/1.2, output.get(4), EPSILON);
    }

    @Test
    public void testMinCoverage()
    {
        final List<Double> input = Arrays.asList(1.0, 1.5, 2.0, -1.0, -1.0);
        maxWindowDistance = 1;
        minWindowCoverage = 3;
        List<Double> output = doNormalisation(input, 1.5);
        assertEquals(input.size(), output.size());
        assertEquals(1.0, output.get(0), EPSILON);
        assertEquals(1.5, output.get(1), EPSILON);
        assertEquals(2.0, output.get(2), EPSILON);
        assertEquals(-1.0, output.get(3), EPSILON);
        assertEquals(-1.0, output.get(4), EPSILON);
    }

    private List<Double> doNormalisation(List<Double> inputs, double expectedRatio)
    {
        int position = 1;
        List<BamRatio> ratios = new ArrayList<>(inputs.size());
        for(Double input : inputs)
        {
            ratios.add(br(position, input));
        }
        DiploidRatioNormaliser normaliser = new DiploidRatioNormaliser(maxWindowDistance, minWindowCoverage);
        ratios.forEach(normaliser::recordRatio);
        normaliser.dataCollectionFinished();
        normaliser.setExpectedRatio(expectedRatio);

        List<Double> results = new ArrayList<>();
        ratios.forEach(b -> results.add(normaliser.normalise(b)));
        return results;
    }

    private BamRatio br(int pos, double ratio)
    {
        return br(_1, pos, ratio, 0.45, true);
    }

    private static void assertRatio(final double input, final double output, double median)
    {
        assertEquals(input / median, output, EPSILON);
    }
}
