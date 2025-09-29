package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.cobalt.count.DepthReading;

import org.junit.Test;

public class BamRatioTest
{
    DepthReading readDepth = new DepthReading("1", 1001, 82, 0.49);

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
    }

    @Test
    public void gcNormaliseNegativeValue()
    {
        BamRatio ratio = new BamRatio(_1, readDepth, true);
        ratio.normaliseForGc(-1.0);
        checkBlanked(ratio);
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
    public void toStringTest()
    {
        assertTrue(new BamRatio(_1, readDepth, true).toString().contains("82"));
    }

    private void checkBlanked(BamRatio ratio)
    {
        assertEquals(-1.0, ratio.ratio(), 0.001);
        assertEquals(-1.0, ratio.readDepth(), 0.001);
        assertEquals(-1.0, ratio.gcContent(), 0.001);
    }
}
