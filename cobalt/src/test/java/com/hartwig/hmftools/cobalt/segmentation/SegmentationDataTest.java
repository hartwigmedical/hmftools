package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class SegmentationDataTest extends SegmentationTestBase
{
    List<CobaltRatio> ratios = new ArrayList<>();
    Function<CobaltRatio, Double> valueExtractor = cobaltRatio -> 4.0;
    SegmentationData sd;

    @Before
    public void setup()
    {
        ratios = new ArrayList<>();
        ratios.add(ratioForValue(_1, 1, 0.5));
        ratios.add(ratioForValue(_1, 1001, 0.7));
        ratios.add(ratioForValue(_1, 2001, 0.6));
        sd = new SegmentationData(ratios, valueExtractor);
    }

    @After
    public void cleanup()
    {
    }

    @Test
    public void writeSegmentationFileTest()
    {
        assertEquals(3, sd.count());
    }

    @Test
    public void emptyTest()
    {
        assertFalse(sd.isEmpty());
        sd = new SegmentationData(new ArrayList<>(), valueExtractor);
        assertTrue(sd.isEmpty());
    }

    @Test
    public void base2LogsOfValuesAreUsed()
    {
        assertEquals(3, sd.valuesForSegmentation().length);
        assertEquals(2.0, sd.valuesForSegmentation()[0], 0.001);
        assertEquals(2.0, sd.valuesForSegmentation()[1], 0.001);
        assertEquals(2.0, sd.valuesForSegmentation()[2], 0.001);
    }

    @Test
    public void rawValues()
    {
        ratios = new ArrayList<>();
        ratios.add(ratio(_1, 1, 10.0));
        ratios.add(ratio(_1, 1001, 20.0));
        sd = new SegmentationData(ratios, CobaltRatio::tumorGCRatio);
        assertEquals(2, sd.rawValues().length);
        assertEquals(10.0, sd.rawValues()[0], 0.001);
        assertEquals(20.0, sd.rawValues()[1], 0.001);
    }

    @Test
    public void lowerBoundOnLoggedValues()
    {
        ratios = new ArrayList<>();
        ratios.add(ratio(_1, 1, 0.0009999));
        ratios.add(ratio(_1, 1001, 0.001));
        ratios.add(ratio(_1, 2001, 0.01));
        sd = new SegmentationData(ratios, CobaltRatio::tumorGCRatio);
        assertEquals(3, sd.valuesForSegmentation().length);
        assertEquals(-9.966, sd.valuesForSegmentation()[0], 0.01);
        assertEquals(log2(0.001), sd.valuesForSegmentation()[1], 0.001);
        assertEquals(log2(0.01), sd.valuesForSegmentation()[2], 0.001);
    }

    private double log2(double value)
    {
        return Math.log(value) / Math.log(2);
    }
}
