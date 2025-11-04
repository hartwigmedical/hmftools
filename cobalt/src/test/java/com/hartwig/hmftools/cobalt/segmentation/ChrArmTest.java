package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ChrArmTest
{
    @Test
    public void compareTest()
    {
        assertEquals(0, new ChrArm(_1, Arm.P).compareTo(new ChrArm(_1, Arm.P)));
        assertTrue(new ChrArm(_1, Arm.P).compareTo(new ChrArm(_1, Arm.Q)) < 0);
        assertTrue(new ChrArm(_2, Arm.P).compareTo(new ChrArm(_1, Arm.Q)) > 0);
        assertTrue(new ChrArm(_2, Arm.P).compareTo(new ChrArm(_1, Arm.P)) > 0);
    }
}
