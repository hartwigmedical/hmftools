package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.junit.Assert;
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

    @Test
    public void comparisonTest()
    {
        ArmEvidence ae1p = new ArmEvidence(new ChrArm(_1, Arm.P));
        ArmEvidence ae1q = new ArmEvidence(new ChrArm(_1, Arm.Q));
        ArmEvidence ae2p = new ArmEvidence(new ChrArm(_2, Arm.P));
        ArmEvidence ae2q = new ArmEvidence(new ChrArm(_2, Arm.Q));
        ArmEvidence ae3p = new ArmEvidence(new ChrArm(_3, Arm.P));
        ArmEvidence ae3q = new ArmEvidence(new ChrArm(_3, Arm.Q));

        setValues(ae1p, 1, 10);
        setValues(ae2q, 10, 100);
        setValues(ae1q, 10, 20);
        setValues(ae2p, 20, 20);
        setValues(ae3p, 5, 200);
        setValues(ae3q, 5, 200);

        Assert.assertTrue(ae1p.compareTo(ae2p) < 0);
        Assert.assertTrue(ae1p.compareTo(ae2q) < 0);
        Assert.assertTrue(ae1p.compareTo(ae1q) < 0);
        Assert.assertTrue(ae1q.compareTo(ae2p) < 0);
        Assert.assertTrue(ae3p.compareTo(ae3q) < 0);
    }

    static void setValues(ArmEvidence evidence, int hits, int totalPoints)
    {
        Preconditions.checkArgument(hits <= totalPoints);

        for(int i = 0; i < hits; i++)
        {
            evidence.register(true);
        }
        for(int i = hits; i < totalPoints; i++)
        {
            evidence.register(false);
        }
    }
}
