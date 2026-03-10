package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.purity.CategoryEvidenceTest.setValues;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.junit.Test;

public class CategoryEvidenceIntegralTest
{
    @Test
    public void totalHitsTest()
    {
        CategoryEvidence<ChrArm> ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);
        CategoryEvidence<ChrArm> ae2p = ae(HumanChromosome._2, Arm.P, 5, 10);
        CategoryEvidence<ChrArm> ae3p = ae(HumanChromosome._3, Arm.P, 2, 20);

        CategoryEvidenceIntegral<ChrArm> integral = new CategoryEvidenceIntegral<>(Set.of(ae1p, ae2p, ae3p));
        assertEquals(17, integral.totalHits());
    }

    @Test
    public void totalPointsTest()
    {
        CategoryEvidence<ChrArm> ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);
        CategoryEvidence<ChrArm> ae2p = ae(HumanChromosome._2, Arm.P, 5, 10);
        CategoryEvidence<ChrArm> ae3p = ae(HumanChromosome._3, Arm.P, 2, 20);

        CategoryEvidenceIntegral<ChrArm> integral = new CategoryEvidenceIntegral<>(Set.of(ae1p, ae2p, ae3p));
        assertEquals(50, integral.totalPoints());
    }

    @Test
    public void integralWithOneArmTest()
    {
        CategoryEvidence<ChrArm> ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);

        CategoryEvidenceIntegral<ChrArm> integral = new CategoryEvidenceIntegral<>(Set.of(ae1p));
        assertEquals(100.0, integral.value(), 0.001);
    }

    @Test
    public void integralWithThreeArmsTest()
    {
        CategoryEvidence<ChrArm> ae1p = ae(HumanChromosome._1, Arm.P, 10, 20);
        CategoryEvidence<ChrArm> ae2p = ae(HumanChromosome._2, Arm.P, 5, 10);
        CategoryEvidence<ChrArm> ae3p = ae(HumanChromosome._3, Arm.P, 2, 20);
        // 2 * 20 * 0.5 + (30 * 2) + (30 * 15 * 0.5)
        CategoryEvidenceIntegral<ChrArm> integral = new CategoryEvidenceIntegral<>(Set.of(ae1p, ae2p, ae3p));
        assertEquals(305.0, integral.value(), 0.001);
    }

    @Test
    public void integralWithMultipleArmsTest()
    {
        CategoryEvidence<ChrArm> ae1p = ae(HumanChromosome._1, Arm.P, 2, 10);  // 0.2
        CategoryEvidence<ChrArm> ae1q = ae(HumanChromosome._1, Arm.Q, 3, 10); // 0.3
        CategoryEvidence<ChrArm> ae2p = ae(HumanChromosome._2, Arm.P, 1, 8);  // 0.125
        CategoryEvidence<ChrArm> ae2q = ae(HumanChromosome._2, Arm.Q, 2, 12); // ~0.17
        CategoryEvidence<ChrArm> ae3p = ae(HumanChromosome._3, Arm.P, 1, 10); // 0.1

        // (1*10*0.5) + (1*8 + 1*8*0.5) + (2*12 + 2*12*0.5) + (4*10 + 2*10*0.5) + (6*10 + 3*10*0.5)
        // = 5 + (8 + 4) + (24 + 12) + (40 + 10) + (60 + 15) = 182
        CategoryEvidenceIntegral<ChrArm> integral = new CategoryEvidenceIntegral<>(Set.of(ae1p, ae1q, ae2p, ae2q, ae3p));
        assertEquals(178.0, integral.value(), 0.001);
    }

    private CategoryEvidence<ChrArm> ae(HumanChromosome chromosome, Arm arm, int hits, int pointsInArm)
    {
        CategoryEvidence<ChrArm> evidence = new CategoryEvidence<>(new ChrArm(chromosome, arm));
        setValues(evidence, hits, pointsInArm);
        return evidence;
    }
}
