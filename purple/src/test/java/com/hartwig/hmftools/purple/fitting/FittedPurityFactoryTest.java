package com.hartwig.hmftools.purple.fitting;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.purple.fitting.FittedPurityFactory;

import org.junit.Test;

public class FittedPurityFactoryTest
{
    private static final double EPSILON = 1e-10;

    @Test
    public void testPloidyRange()
    {
        final List<Double> oneToThree = FittedPurityFactory.ploidyRange(1, 3);
        assertEquals(101, oneToThree.size());
        assertEquals(1d, oneToThree.get(0), EPSILON);
        assertEquals(3d, oneToThree.get(100), EPSILON);

        final List<Double> threeToFive = FittedPurityFactory.ploidyRange(3, 5);
        assertEquals(41, threeToFive.size());
        assertEquals(3d, threeToFive.get(0), EPSILON);
        assertEquals(5d, threeToFive.get(40), EPSILON);

        final List<Double> fiveOnwards = FittedPurityFactory.ploidyRange(5, 8);
        assertEquals(31, fiveOnwards.size());
        assertEquals(5d, fiveOnwards.get(0), EPSILON);
        assertEquals(8d, fiveOnwards.get(30), EPSILON);

        final List<Double> all = FittedPurityFactory.ploidyRange(1, 8);
        assertEquals(171, all.size());
        assertEquals(1d, all.get(0), EPSILON);
        assertEquals(3d, all.get(100), EPSILON);
        assertEquals(5d, all.get(140), EPSILON);
        assertEquals(8d, all.get(170), EPSILON);
    }

    @Test
    public void testFixedPloidyRange()
    {
        double fixedPloidy = 0.3456;
        final List<Double> fixed = FittedPurityFactory.ploidyRange(fixedPloidy, fixedPloidy);
        assertEquals(1, fixed.size());
        assertEquals(fixedPloidy, fixed.get(0), EPSILON);
    }
}
