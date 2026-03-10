package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.junit.Test;

public class CandidatePeakTest extends PurityTestBase
{

    @Test
    public void shouldGoInBucketToWhichItIsClosest()
    {
        CandidatePeak peak = new CandidatePeak(0.006, 0.01);
        peak.test(evidenceWithDepthAndAltCount(_1, 1000, 1746, 7));
        assertEquals(1, peak.heterozygousEvidencePoints().size());
        assertEquals(0, peak.homozygousEvidencePoints().size());

        peak = new CandidatePeak(0.006, 0.01);
        peak.test(evidenceWithDepthAndAltCount(_1, 1000, 1746, 8));
        assertEquals(0, peak.heterozygousEvidencePoints().size());
        assertEquals(1, peak.homozygousEvidencePoints().size());
    }

    @Test
    public void hasSufficientDepthForEventDetectionTest()
    {
        // 16th percentile of het peak @ 2 AD
        CandidatePeak level35 = new CandidatePeak(0.35);
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(0)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(1)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(20)));
        assertFalse(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(25)));
        assertTrue(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(26)));
        assertTrue(level35.hasSufficientDepthForEventDetection(evidenceWithDepth(30)));

        CandidatePeak level10 = new CandidatePeak(0.1);
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
    public void altAndRefDepthAreUsedToCalculateVafTest()
    {
        CandidatePeak candidatePeak = new CandidatePeak(0.2);
        PositionEvidence pe = new PositionEvidence("1", 1000, "A", "T");
        pe.AltSupport = 10;
        pe.RefSupport = 40;
        pe.ReadDepth = 10000;
        candidatePeak.test(pe);
        assertTrue(candidatePeak.allCapturedPoints().contains(pe));
    }

    @Test
    public void altAndRefDepthAreEffectiveReadDepthTest()
    {
        CandidatePeak candidatePeak = new CandidatePeak(0.02);
        PositionEvidence pe = new PositionEvidence("1", 1000, "A", "T");
        pe.AltSupport = 10;
        pe.RefSupport = 40;
        pe.ReadDepth = 10000;
        assertFalse(candidatePeak.hasSufficientDepthForEventDetection(pe));
    }

    @Test
    public void numberOfEventsCapturedTest()
    {
        CandidatePeak level = new CandidatePeak(0.1);
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
        CandidatePeak level = new CandidatePeak(0.1);
        assertEquals(0.1, level.vaf(), 0.0001);
    }

    @Test
    public void homozygousProportionTest()
    {
        CandidatePeak level = new CandidatePeak(0.1);
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
        CandidatePeak level = new CandidatePeak(0.1, 0.02);
        level.test(evidenceWithDepthAndAltCount(_3, 1_000_000, 1000, 100));
        level.test(evidenceWithDepthAndAltCount(_3, 2_000_000, 1000, 50));
        level.test(evidenceWithDepthAndAltCount(_3, 2_001_000, 1000, 50));
        level.test(evidenceWithDepthAndAltCount(_3, 3_000_000, 1000, 10));
        assertEquals("VafLevel{vaf=0.10, step=0.02, tested: 4, homozygous: 1, heterozygous: 2}", level.toString());
    }

    private void checkNotCaptured(double vafLevel, int readDepth, int altDepth)
    {
        CandidatePeak level = new CandidatePeak(vafLevel);
        checkNotCaptured(readDepth, altDepth, level);
    }

    private void checkNotCaptured(double vafLevel, double gap, int readDepth, int altDepth)
    {
        CandidatePeak level = new CandidatePeak(vafLevel, gap);
        checkNotCaptured(readDepth, altDepth, level);
    }

    private void checkNotCaptured(final int readDepth, final int altDepth, final CandidatePeak level)
    {
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.homozygousEvidencePoints().isEmpty());
        assertTrue(level.heterozygousEvidencePoints().isEmpty());
    }

    private void checkCapturedHom(double vafLevel, int readDepth, int altDepth)
    {
        CandidatePeak level = new CandidatePeak(vafLevel);
        checkCapturedHom(readDepth, altDepth, level);
    }

    private void checkCapturedHom(double vafLevel, double gap, int readDepth, int altDepth)
    {
        CandidatePeak level = new CandidatePeak(vafLevel, gap);
        checkCapturedHom(readDepth, altDepth, level);
    }

    private void checkCapturedHom(final int readDepth, final int altDepth, final CandidatePeak level)
    {
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.homozygousEvidencePoints().contains(positionEvidence));
    }

    private void checkCapturedHet(double vafLevel, int readDepth, int altDepth)
    {
        CandidatePeak level = new CandidatePeak(vafLevel);
        checkCapturedHet(readDepth, altDepth, level);
    }

    private void checkCapturedHet(double vafLevel, double gap, int readDepth, int altDepth)
    {
        CandidatePeak level = new CandidatePeak(vafLevel, gap);
        checkCapturedHet(readDepth, altDepth, level);
    }

    private void checkCapturedHet(final int readDepth, final int altDepth, final CandidatePeak level)
    {
        final PositionEvidence positionEvidence = evidenceWithDepthAndAltCount(readDepth, altDepth);
        level.test(positionEvidence);
        assertTrue(level.heterozygousEvidencePoints().contains(positionEvidence));
    }
}
