package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.junit.Test;

public class ArmEvidenceTest
{
    @Test
    public void equalsTest()
    {
        ArmEvidence evidence1 = new ArmEvidence(new ChrArm(_1, Arm.P));
        ArmEvidence evidence2 = new ArmEvidence(new ChrArm(_1, Arm.P));
        ArmEvidence evidence3 = new ArmEvidence(new ChrArm(_1, Arm.Q));
        ArmEvidence evidence4 = new ArmEvidence(new ChrArm(_2, Arm.P));

        assertEquals(evidence1, evidence2);
        assertNotEquals(evidence1, evidence3);
        assertNotEquals(evidence1, evidence4);
    }

    @Test
    public void hashCodeTest()
    {
        ArmEvidence evidence1 = new ArmEvidence(new ChrArm(_1, Arm.P));
        ArmEvidence evidence2 = new ArmEvidence(new ChrArm(_1, Arm.P));

        assertEquals(evidence1.hashCode(), evidence2.hashCode());
    }

    @Test
    public void ratioTest()
    {
        ArmEvidence evidence = new ArmEvidence(new ChrArm(_1, Arm.P));
        assertEquals(Double.NaN, evidence.ratio(), 0.0001);
        evidence.register(false);
        assertEquals(0.0, evidence.ratio(), 0.001);
        evidence.register(true);
        assertEquals(0.5, evidence.ratio(), 0.001);
    }
}
