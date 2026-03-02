package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.amber.AmberBase.A;
import static com.hartwig.hmftools.common.amber.AmberBase.C;
import static com.hartwig.hmftools.common.amber.AmberBase.G;
import static com.hartwig.hmftools.common.amber.AmberBase.T;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.junit.Test;

public class VafLevelTest extends PurityTestBase
{
    private static final int CENTRE = 10_000_000;
    private final ChrArmLocator Locator = position ->
    {
        Arm arm = position.position() < CENTRE ? Arm.P : Arm.Q;
        return new ChrArm(position.chr(), arm);
    };

    @Test
    public void hasSufficientDepthForEventDetectionTest()
    {
        BinomialDistribution binomial = new BinomialDistribution(34, 0.35);
        var v = binomial.cumulativeProbability(2);
        // 16th percentile of het peak @ 2 AD
        VafLevel level35 = new VafLevel(0.35);
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(0)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(11)));
        assertTrue(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(12)));

        VafLevel level10 = new VafLevel(0.1);
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(0)));
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(10)));
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(44)));
        assertTrue(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(45)));
        assertTrue(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(125)));

        VafLevel level1 = new VafLevel(0.01);
        assertFalse(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(0)));
        assertFalse(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(100)));
        assertFalse(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(450)));
        assertFalse(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(461)));
        assertTrue(level1.hasSufficientDepthForEventDetection(evidenceWithDepth(462)));
    }

    @Test
    public void testTest()
    {
        checkNotCaptured(0.1, 1000, 20);
        checkNotCaptured(0.1, 1000, 44);
        checkCapturedHet(0.1, 1000, 45);
        checkCapturedHet(0.1, 1000, 54);
        checkNotCaptured(0.1, 1000, 55);
        checkNotCaptured(0.1, 1000, 80);
        //        checkNotCaptured(0.1, 1000, 94); TODO
        checkCapturedHom(0.1, 1000, 95);
        checkCapturedHom(0.1, 1000, 100);
        checkCapturedHom(0.1, 1000, 105);
        checkNotCaptured(0.1, 1000, 106);

        checkNotCaptured(0.30, 100, 20);
        checkNotCaptured(0.30, 100, 26);
        checkCapturedHom(0.30, 100, 27);
        checkCapturedHom(0.30, 100, 30);
        checkCapturedHom(0.30, 100, 32);
        checkNotCaptured(0.30, 100, 34);
    }

    @Test
    public void evenCaptureAcrossChromosomeArmsTest()
    {
        VafLevel level = new VafLevel(0.1);
        // 10 reads per chromosome arm, 2 of which are captured.
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            int startPos_P = 1_000_000;
            int startPos_Q = CENTRE + 1_000_000;
            for(int i = 0; i < 8; i++)
            {
                level.test(evidenceWithDepthAndAltCount(chromosome, startPos_P + i * 1000, 1000, 10));
                level.test(evidenceWithDepthAndAltCount(chromosome, startPos_Q + i * 1000, 1000, 10));
            }
            level.test(evidenceWithDepthAndAltCount(chromosome, startPos_P + 8 * 1000, 1000, 100));
            level.test(evidenceWithDepthAndAltCount(chromosome, startPos_Q + 8 * 1000, 1000, 100));
            level.test(evidenceWithDepthAndAltCount(chromosome, startPos_P + 9 * 1000, 1000, 100));
            level.test(evidenceWithDepthAndAltCount(chromosome, startPos_Q + 9 * 1000, 1000, 100));
        }
        assertEquals(1.0, level.perArmConsistencyFactor(Locator), 0.0001);
    }

    @Test
    public void allCapturedEventsAreInASingleChromosomeArmTest()
    {
        VafLevel level = new VafLevel(0.1);
        // 8 reads per chromosome arm, none of which are captured.
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            int startPos_P = 1_000_000;
            int startPos_Q = CENTRE + 1_000_000;
            for(int i = 0; i < 8; i++)
            {
                level.test(evidenceWithDepthAndAltCount(chromosome, startPos_P + i * 1000, 1000, 10));
                level.test(evidenceWithDepthAndAltCount(chromosome, startPos_Q + i * 1000, 1000, 10));
            }
        }
        // Two events in 1P that are captured
        level.test(evidenceWithDepthAndAltCount(_1, 2_000_000, 1000, 100));
        level.test(evidenceWithDepthAndAltCount(_1, 2_001_000, 1000, 100));

        // 10/(46 * 8 + 2)
        assertEquals(0.0259, level.perArmConsistencyFactor(Locator), 0.0001);
    }

    @Test
    public void eventsEvenlyCapturedAccordingToMutationTypeTest()
    {
        VafLevel level = new VafLevel(0.2);
        int position = 1_000_000;
        // For each possible SNV create 10 reads, of which 1 is captured.
        for(AmberBase ref : AmberBase.values())
        {
            for(AmberBase alt : AmberBase.values())
            {
                if(ref == alt)
                {
                    continue;
                }
                for(int i = 0; i < 9; i++)
                {
                    position += 1000;
                    level.test(evidenceWithDepthAndAltCount(ref, alt, position, 1000, 10));
                }
                position += 1000;
                level.test(evidenceWithDepthAndAltCount(ref, alt, position, 1000, 200));
            }
        }
        // Each of the 6 canonical SNV types has 10 reads, of which 1 is in the peak's hom range.
        assertEquals(1.0, level.perMutationTypeConsistencyFactor(), 0.0001);
    }

    @Test
    public void eventsOnlyCapturedInCAMutationsTest()
    {
        VafLevel level = new VafLevel(0.3);
        int position = 1_000_000;
        for(AmberBase ref : AmberBase.values())
        {
            for(AmberBase alt : AmberBase.values())
            {
                if(ref == alt)
                {
                    continue;
                }
                for(int i = 0; i < 4; i++)
                {
                    position += 1000;
                    level.test(evidenceWithDepthAndAltCount(ref, alt, position, 1000, 10));
                }
                position += 1000;
                if(ref == C && alt == A || ref == G && alt == T)
                {
                    level.test(evidenceWithDepthAndAltCount(ref, alt, position, 1000, 300));
                }
                else
                {
                    level.test(evidenceWithDepthAndAltCount(ref, alt, position, 1000, 10));
                }
            }
        }
        // Each of the 6 canonical SNV types has 10 reads. For C>A, 2 of these are captured in the value's hom peak.
        // So expected = 1/6.
        assertEquals(0.1667, level.perMutationTypeConsistencyFactor(), 0.0001);
    }

    @Test
    public void numberOfEventsCapturedTest()
    {
        VafLevel level = new VafLevel(0.1);
        assertEquals(0, level.numberOfCapturedEvidencePoints());
        level.test(evidenceWithDepthAndAltCount(_3, 1_000_000, 1000, 100));
        assertEquals(1, level.numberOfCapturedEvidencePoints());
        level.test(evidenceWithDepthAndAltCount(_3, 2_000_000, 1000, 50));
        assertEquals(2, level.numberOfCapturedEvidencePoints());
        level.test(evidenceWithDepthAndAltCount(_3, 3_000_000, 1000, 10));
        assertEquals(2, level.numberOfCapturedEvidencePoints());
    }

    @Test
    public void vafTest()
    {
        VafLevel level = new VafLevel(0.1);
        assertEquals(0.1, level.vaf(), 0.0001);
    }

    @Test
    public void homozygousProportionTest()
    {
        VafLevel level = new VafLevel(0.1);
        assertEquals(Double.NaN, level.homozygousProportion(), 0.0001);

        level.test(evidenceWithDepthAndAltCount(_3, 3_000_000, 1000, 10));
        assertEquals(Double.NaN, level.homozygousProportion(), 0.0001);

        level.test(evidenceWithDepthAndAltCount(_3, 4_000_000, 1000, 50)); // het
        assertEquals(0.0, level.homozygousProportion(), 0.0001);

        level.test(evidenceWithDepthAndAltCount(_3, 5_000_000, 1000, 50)); // het
        assertEquals(0.0, level.homozygousProportion(), 0.0001);

        level.test(evidenceWithDepthAndAltCount(_3, 6_000_000, 1000, 100)); // hom
        assertEquals(0.3333, level.homozygousProportion(), 0.0001);

        level.test(evidenceWithDepthAndAltCount(_3, 7_000_000, 1000, 100)); // hom
        assertEquals(0.5, level.homozygousProportion(), 0.0001);
    }

    private void checkNotCaptured(double vafLevel, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel);
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.homozygousEvidencePoints().isEmpty());
        assertTrue(level.heterozygousEvidencePoints().isEmpty());
    }

    private void checkCapturedHom(double vafLevel, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel);
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.homozygousEvidencePoints().contains(positionEvidence));
    }

    private void checkCapturedHet(double vafLevel, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel);
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.heterozygousEvidencePoints().contains(positionEvidence));
    }
}
