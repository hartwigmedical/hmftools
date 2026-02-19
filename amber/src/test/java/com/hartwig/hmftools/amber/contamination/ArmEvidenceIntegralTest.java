package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.amber.contamination.ArmEvidenceTest.setValues;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.junit.Test;

public class ArmEvidenceIntegralTest
{
    @Test
    public void totalHitsTest()
    {
        ArmEvidence ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);
        ArmEvidence ae2p = ae(HumanChromosome._2, Arm.P, 5, 10);
        ArmEvidence ae3p = ae(HumanChromosome._3, Arm.P, 2, 20);

        ArmEvidenceIntegral integral = new ArmEvidenceIntegral(Set.of(ae1p, ae2p, ae3p));
        assertEquals(17, integral.totalHits());
    }

    @Test
    public void totalPointsTest()
    {
        ArmEvidence ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);
        ArmEvidence ae2p = ae(HumanChromosome._2, Arm.P, 5, 10);
        ArmEvidence ae3p = ae(HumanChromosome._3, Arm.P, 2, 20);

        ArmEvidenceIntegral integral = new ArmEvidenceIntegral(Set.of(ae1p, ae2p, ae3p));
        assertEquals(50, integral.totalPoints());
    }

    @Test
    public void integralWithOneArmTest()
    {
        ArmEvidence ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);

        ArmEvidenceIntegral integral = new ArmEvidenceIntegral(Set.of(ae1p));
        assertEquals(100.0, integral.value(), 0.001);
    }

    @Test
    public void integralWithThreeArmsTest()
    {
        ArmEvidence ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);
        ArmEvidence ae2p = ae(HumanChromosome._2, Arm.P, 5, 10);
        ArmEvidence ae3p = ae(HumanChromosome._3, Arm.P, 2, 20);
        // 2 * 20 * 0.5 + (30 * 2) + (30 * 15 * 0.5)
        ArmEvidenceIntegral integral = new ArmEvidenceIntegral(Set.of(ae1p, ae2p, ae3p));
        assertEquals(305.0, integral.value(), 0.001);
    }

    @Test
    public void integralWithMultipleArmsTest()
    {
        ArmEvidence ae1p = ae(HumanChromosome._1, Arm.P, 2, 10);  // 0.2
        ArmEvidence ae1q = ae(HumanChromosome._1, Arm.Q, 3, 10); // 0.3
        ArmEvidence ae2p = ae(HumanChromosome._2, Arm.P, 1, 8);  // 0.125
        ArmEvidence ae2q = ae(HumanChromosome._2, Arm.Q, 2, 12); // ~0.17
        ArmEvidence ae3p = ae(HumanChromosome._3, Arm.P, 1, 10); // 0.1

        // (1*10*0.5) + (1*8 + 1*8*0.5) + (2*12 + 2*12*0.5) + (4*10 + 2*10*0.5) + (6*10 + 3*10*0.5)
        // = 5 + (8 + 4) + (24 + 12) + (40 + 10) + (60 + 15) = 182
        ArmEvidenceIntegral integral = new ArmEvidenceIntegral(Set.of(ae1p, ae1q, ae2p, ae2q, ae3p));
        assertEquals(178.0, integral.value(), 0.001);
    }

    private ArmEvidence ae(HumanChromosome chromosome, Arm arm, int hits, int pointsInArm)
    {
        ArmEvidence evidence = new ArmEvidence(new ChrArm(chromosome, arm));
        setValues(evidence, hits, pointsInArm);
        return evidence;
    }
}
