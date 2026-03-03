package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.amber.AmberBase.A;
import static com.hartwig.hmftools.common.amber.AmberBase.C;
import static com.hartwig.hmftools.common.amber.AmberBase.G;
import static com.hartwig.hmftools.common.amber.AmberBase.N;
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
        // 16th percentile of het peak @ 2 AD
        VafLevel level35 = new VafLevel(0.35);
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(0)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(20)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(25)));
        assertTrue(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(26)));
        assertTrue(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(30)));

        VafLevel level10 = new VafLevel(0.1);
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(0)));
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(91)));
        assertTrue(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(92)));
        assertTrue(level10.hasSufficientDepthForEventDetection(evidenceWithDepth(1000)));
    }

    @Test
    public void testCaptureByBinomialTest()
    {
        checkNotCaptured(0.1, 1000, 20);
        checkNotCaptured(0.1, 1000, 42);
        checkCapturedHet(0.1, 1000, 50);
        checkCapturedHet(0.1, 1000, 56);
        checkNotCaptured(0.1, 1000, 57);
        checkNotCaptured(0.1, 1000, 79);
        checkNotCaptured(0.1, 1000, 90);
        checkCapturedHom(0.1, 1000, 91);
        checkCapturedHom(0.1, 1000, 100);
        checkCapturedHom(0.1, 1000, 108);
        checkNotCaptured(0.1, 1000, 110);

        checkNotCaptured(0.30, 100, 5);
        checkNotCaptured(0.30, 100, 10);
        checkCapturedHet(0.30, 100, 11);
        checkCapturedHet(0.30, 100, 15);
        checkCapturedHet(0.30, 100, 18);
        checkNotCaptured(0.30, 100, 19);
        checkNotCaptured(0.30, 100, 20);
        checkNotCaptured(0.30, 100, 24);
        checkCapturedHom(0.30, 100, 25);
        checkCapturedHom(0.30, 100, 30);
        checkCapturedHom(0.30, 100, 34);
        checkNotCaptured(0.30, 100, 35);
    }

    @Test
    public void captureByNextLevelStepTest()
    {
        new BinomialDistribution(1000, 0.1).cumulativeProbability(90);
        checkNotCaptured(0.1, 0.018, 1000, 40);
        checkNotCaptured(0.1, 0.019, 1000, 40);
        checkNotCaptured(0.1, 0.020, 1000, 40);
        checkCapturedHet(0.1, 0.020001, 1000, 40);

        checkNotCaptured(0.1, 0.018, 1000, 60);
        checkNotCaptured(0.1, 0.019, 1000, 60);
        checkNotCaptured(0.1, 0.019999, 1000, 60);
        checkCapturedHet(0.1, 0.02, 1000, 60);

        checkNotCaptured(0.1, 0.02, 1000, 79);
        checkNotCaptured(0.1, 0.019, 1000, 80);
        checkNotCaptured(0.1, 0.02, 1000, 80);
        checkCapturedHom(0.1, 0.020001, 1000, 80);
        checkCapturedHom(0.1, 0.02, 1000, 120);
        checkNotCaptured(0.1, 0.02, 1000, 121);
    }

    @Test
    public void captureByNextLevelStepWhenVafCloseTo1Test()
    {
        new BinomialDistribution(1000, 0.1).cumulativeProbability(90);
        checkNotCaptured(0.1, 0.018, 1000, 960);
        checkNotCaptured(0.1, 0.019, 1000, 960);
        checkNotCaptured(0.1, 0.020, 1000, 960);
        checkCapturedHet(0.1, 0.020001, 1000, 960);

        checkNotCaptured(0.1, 0.018, 1000, 940);
        checkNotCaptured(0.1, 0.019, 1000, 940);
        checkNotCaptured(0.1, 0.019999, 1000, 940);
        checkCapturedHet(0.1, 0.02, 1000, 940);

        checkNotCaptured(0.1, 0.02, 1000, 921);
        checkNotCaptured(0.1, 0.019, 1000, 920);
        checkNotCaptured(0.1, 0.02, 1000, 920);
        checkCapturedHom(0.1, 0.020001, 1000, 920);
        checkCapturedHom(0.1, 0.02, 1000, 880);
        checkNotCaptured(0.1, 0.02, 1000, 879);
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
            if(ref == N)
            {
                continue;
            }
            for(AmberBase alt : AmberBase.values())
            {
                if(alt == N || ref == alt)
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
            if(ref == N)
            {
                continue;
            }
            for(AmberBase alt : AmberBase.values())
            {
                if(alt == N || ref == alt)
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

    @Test
    public void toStringTest()
    {
        VafLevel level = new VafLevel(0.1, 0.02);
        level.test(evidenceWithDepthAndAltCount(_3, 1_000_000, 1000, 100));
        level.test(evidenceWithDepthAndAltCount(_3, 2_000_000, 1000, 50));
        level.test(evidenceWithDepthAndAltCount(_3, 2_001_000, 1000, 50));
        level.test(evidenceWithDepthAndAltCount(_3, 3_000_000, 1000, 10));
        assertEquals("VafLevel{vaf=0.10, step=0.02, tested: 4, homozygous: 1, heterozygous: 2}", level.toString());
    }

    private void checkNotCaptured(double vafLevel, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel);
        checkNotCaptured(readDepth, altDepth, level);
    }

    private void checkNotCaptured(double vafLevel, double gap, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel, gap);
        checkNotCaptured(readDepth, altDepth, level);
    }

    private void checkNotCaptured(final int readDepth, final int altDepth, final VafLevel level)
    {
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.homozygousEvidencePoints().isEmpty());
        assertTrue(level.heterozygousEvidencePoints().isEmpty());
    }

    private void checkCapturedHom(double vafLevel, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel);
        checkCapturedHom(readDepth, altDepth, level);
    }

    private void checkCapturedHom(double vafLevel, double gap, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel, gap);
        checkCapturedHom(readDepth, altDepth, level);
    }

    private void checkCapturedHom(final int readDepth, final int altDepth, final VafLevel level)
    {
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.homozygousEvidencePoints().contains(positionEvidence));
    }

    private void checkCapturedHet(double vafLevel, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel);
        checkCapturedHet(readDepth, altDepth, level);
    }

    private void checkCapturedHet(double vafLevel, double gap, int readDepth, int altDepth)
    {
        VafLevel level = new VafLevel(vafLevel, gap);
        checkCapturedHet(readDepth, altDepth, level);
    }

    private void checkCapturedHet(final int readDepth, final int altDepth, final VafLevel level)
    {
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.heterozygousEvidencePoints().contains(positionEvidence));
    }
}
