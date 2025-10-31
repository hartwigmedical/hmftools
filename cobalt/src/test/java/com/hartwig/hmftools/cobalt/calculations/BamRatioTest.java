package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.cobalt.count.DepthReading;

import org.junit.Test;

public class BamRatioTest
{
    DepthReading readDepth = new DepthReading("1", 1001, 82, 0.49);

    @Test
    public void depth0Test()
    {
        DepthReading depth0Reading = new DepthReading("1", 1001, 0.0, Double.NaN);
        BamRatio br = new BamRatio(_1, depth0Reading, true);

        assertEquals(0.0, br.readDepth(), 0.0001);
        assertEquals(0.0, br.ratio(), 0.0001);
        assertEquals(0.0, br.gcContent(), 0.0001);
        assertEquals(-1.0, br.getDiploidAdjustedRatio(), 0.0001);

        // Depth 0 readings are invariant under normalisation and enrichment.
        br.normaliseByMean(3.4);
        assertEquals(0.0, br.readDepth(), 0.0001);
        assertEquals(0.0, br.ratio(), 0.0001);
        assertEquals(0.0, br.gcContent(), 0.0001);
        assertEquals(-1.0, br.getDiploidAdjustedRatio(), 0.0001);

        br.normaliseForGc(45.89);
        assertEquals(0.0, br.readDepth(), 0.0001);
        assertEquals(0.0, br.ratio(), 0.0001);
        assertEquals(0.0, br.gcContent(), 0.0001);
        assertEquals(-1.0, br.getDiploidAdjustedRatio(), 0.0001);

        br.applyEnrichment(34.56);
        assertEquals(0.0, br.readDepth(), 0.0001);
        assertEquals(0.0, br.ratio(), 0.0001);
        assertEquals(0.0, br.gcContent(), 0.0001);
        assertEquals(-1.0, br.getDiploidAdjustedRatio(), 0.0001);
    }

    @Test
    public void constructors()
    {
        BamRatio br = new BamRatio(_Y, 19_001, 123.4, 0.55);
        assertEquals(_Y, br.mChromosome);
        assertEquals(19_001, br.position());
        assertEquals(123.4, br.ratio(), 0.001);
        assertEquals(123.4, br.readDepth(), 0.001);
        assertEquals(0.55, br.gcContent(), 0.001);
        assertEquals(-1.0, br.getDiploidAdjustedRatio(), 0.001);

        br = new BamRatio(_X, 19_001, 123.4, 12.34,0.34);
        assertEquals(_X, br.mChromosome);
        assertEquals(19_001, br.position());
        assertEquals(12.34, br.ratio(), 0.001);
        assertEquals(123.4, br.readDepth(), 0.001);
        assertEquals(0.34, br.gcContent(), 0.001);
        assertEquals(-1.0, br.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void overrideRatio()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(82.0, ratio.ratio(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
        ratio.overrideRatio(1.8);
        assertEquals(1.8, ratio.ratio(), 0.001);
    }

    @Test
    public void overrideRatioSetsIncludedStatus()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, false);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.ratio(), 0.001); // sanity
        ratio.overrideRatio(1.8);
        assertEquals(1.8, ratio.ratio(), 0.001);
    }

    @Test
    public void overrideWithNegativeRatioSetsIncludedStatus()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, false);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.ratio(), 0.001);
        ratio.overrideRatio(-1.0);
        assertEquals(-1.0, ratio.ratio(), 0.001);
    }

    @Test
    public void normaliseByMean()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.normaliseByMean(0.5);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(164, ratio.ratio(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void normalise0ByMean()
    {
        DepthReading rd0 = new DepthReading("1", 1001, 0.0, 0.49);
        BamRatio ratio = new BamRatio(_1, rd0, true);
        assertEquals(0.0, ratio.gcContent(), 0.001);
        ratio.normaliseByMean(0.5);
        assertEquals(0.0, ratio.readDepth(), 0.001);
        assertEquals(0.0, ratio.ratio(), 0.001);
        assertEquals(0.0, ratio.gcContent(), 0.001);
    }

    @Test
    public void normaliseNaNByMean()
    {
        DepthReading rdNaN = new DepthReading("1", 1001, Double.NaN, 0.49);
        BamRatio ratio = new BamRatio(_1, rdNaN, true);
        ratio.normaliseByMean(0.5);
        assertEquals(_1, ratio.mChromosome);
        assertEquals(1001, ratio.Position);
        assertEquals(-1.0, ratio.ratio(), 0.001);
        assertEquals(-1.0, ratio.readDepth(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void normaliseDiploidRatioWhenNotSet()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        assertEquals(-1.0, ratio.getDiploidAdjustedRatio(), 0.001);
        ratio.normaliseDiploidAdjustedRatio(0.5);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(82.0, ratio.ratio(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
        assertEquals(-1.0, ratio.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void normaliseDiploidRatio()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.setDiploidAdjustedRatio(1.8);
        assertEquals(1.8, ratio.getDiploidAdjustedRatio(), 0.001);
        ratio.normaliseDiploidAdjustedRatio(2.0);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(0.9, ratio.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void normaliseDiploidRatioByNaN()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.setDiploidAdjustedRatio(1.8);
        assertEquals(1.8, ratio.getDiploidAdjustedRatio(), 0.001);
        ratio.normaliseDiploidAdjustedRatio(Double.NaN);
        assertEquals(-1.0, ratio.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void normaliseDiploidRatioByZero()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.setDiploidAdjustedRatio(1.8);
        assertEquals(1.8, ratio.getDiploidAdjustedRatio(), 0.001);
        ratio.normaliseDiploidAdjustedRatio(0.0);
        assertEquals(-1.0, ratio.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void diploidRatioWhenReadRatioIsZero()
    {
        BamRatio ratio = new BamRatio(_1, 1001, 0.0, 0.0);
        assertEquals(-1.0, ratio.getDiploidAdjustedRatio(), 0.001);
        ratio.setDiploidAdjustedRatio(1.8);
        assertEquals(0.0, ratio.getDiploidAdjustedRatio(), 0.001);
        ratio.normaliseDiploidAdjustedRatio(10.0);
        assertEquals(0.0, ratio.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void handleNaNDepthInConstructor()
    {
        DepthReading rdNaN = new DepthReading("1", 1001, Double.NaN, 0.49);
        BamRatio ratio = new BamRatio(_1, rdNaN, true);
        assertEquals(_1, ratio.mChromosome);
        assertEquals(1001, ratio.Position);
        assertEquals(-1.0, ratio.ratio(), 0.001);
        assertEquals(-1.0, ratio.readDepth(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void handleNaNGcInConstructor()
    {
        DepthReading rdNaN = new DepthReading("1", 1001, 34, Double.NaN);
        BamRatio ratio = new BamRatio(_1, rdNaN, true);
        assertEquals(_1, ratio.mChromosome);
        assertEquals(1001, ratio.Position);
        assertEquals(34, ratio.ratio(), 0.001);
        assertEquals(34, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.gcContent(), 0.001);
        assertEquals(-1.0, ratio.getDiploidAdjustedRatio(), 0.001);
    }

    @Test
    public void inTargetRegion()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        assertEquals(_1, ratio.mChromosome);
        assertEquals(1001, ratio.Position);
        assertEquals(82, ratio.ratio(), 0.001);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void outsideTargetRegion()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, false);
        assertEquals(_1, ratio.mChromosome);
        assertEquals(1001, ratio.Position);
        checkBlanked(ratio);
    }

    @Test
    public void gcNormalise()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.normaliseForGc(100.0);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(0.82, ratio.ratio(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void gcNormaliseByZero()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.normaliseForGc(0.0);
        checkBlanked(ratio);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void gcNormaliseNegativeValue()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.normaliseForGc(-1.0);
        checkBlanked(ratio);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void applyEnrichment()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.applyEnrichment(0.5);
        assertEquals(82, ratio.readDepth(), 0.001);
        assertEquals(164, ratio.ratio(), 0.001);
        assertEquals(0.49, ratio.gcContent(), 0.001);
    }

    @Test
    public void handleNaNEnrichment()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.applyEnrichment(Double.NaN);
        checkBlanked(ratio);
    }

    @Test
    public void handleNaNEnrichmentForDepth0Window()
    {
        BamRatio ratio = new BamRatio(_1, 1001, 0.0, 0.0);
        ratio.applyEnrichment(Double.NaN);
        checkBlanked(ratio);
    }

    @Test
    public void toStringTest()
    {
        assertTrue(new BamRatio(_1, readDepth, true).toString().contains("82"));
    }

    @Test
    public void position()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        assertEquals(1001, ratio.position());
    }

    private void checkBlanked(BamRatio ratio)
    {
        assertEquals(-1.0, ratio.ratio(), 0.001);
    }

    @Test
    public void testEquals()
    {
        // Same objects should be equal
        BamRatio ratio1 = new BamRatio(_1, readDepth, true);
        BamRatio ratio2 = new BamRatio(_1, readDepth, true);
        assertEquals(ratio1, ratio2);

        // Different chromosome
        BamRatio differentChromosome = new BamRatio(_2, new DepthReading("2", 1001, 82, 0.49), true);
        assertNotEquals(ratio1, differentChromosome);

        // Different position
        BamRatio differentPosition = new BamRatio(_1, new DepthReading("1", 2001, 82, 0.49), true);
        assertNotEquals(ratio1, differentPosition);

        // Different read depth
        BamRatio differentReadDepth = new BamRatio(_1, new DepthReading("1", 1001, 100, 0.49), true);
        assertEquals(ratio1, differentReadDepth);

        // Different GC content
        BamRatio differentGcContent = new BamRatio(_1, new DepthReading("1", 1001, 82, 0.6), true);
        assertEquals(ratio1, differentGcContent);

        // Different included status
        BamRatio differentIncluded = new BamRatio(_1, readDepth, false);
        assertEquals(ratio1, differentIncluded);

        // Different diploid adjusted ratio
        BamRatio differentDiploidRatio = new BamRatio(_1, readDepth, true);
        differentDiploidRatio.setDiploidAdjustedRatio(2.0);
        assertEquals(ratio1, differentDiploidRatio);

        // Null comparison
        assertFalse(ratio1.equals(null));

        // Different class comparison
        assertFalse(ratio1.equals("whatever"));
    }

    @Test
    public void testHashCode()
    {
        BamRatio ratio1 = new BamRatio(_1, readDepth, true);
        BamRatio ratio2 = new BamRatio(_1, readDepth, true);
        assertEquals(ratio1.hashCode(), ratio2.hashCode());
    }
}
