package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.segmentation.PerArmSegmenter;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class TumorRatioSegmenterTest extends SegmentationTestBase
{
    @Before
    public void setup()
    {
        super.setup();
    }

    @After
    public void cleanup()
    {
        super.cleanup();
    }

    @Test
    public void tumorSegmenterTest()
    {
        ratios.put(_1, ratioForValue(_1, 1001, 1.0));
        ratios.put(_1, ratioForValue(_1, 2001, 1.0));
        PerArmSegmenter<CobaltRatio> segmenter = new TumorRatioSegmenter(ratios, locator, 50.0);
        CobaltRatio ratio = new CobaltRatio("1", 1, 2, 2.1, 2.2, 2.3, 3, 3.1, 3.2);
        assertEquals(ratio.tumorGCRatio(), segmenter.value(ratio), 0.00001);
    }

    @Test
    public void isWindowed()
    {
        ratios.put(_1, ratioForValue(_1, 1001, 1.0));
        PerArmSegmenter<CobaltRatio> segmenter = new TumorRatioSegmenter(ratios, locator, 50.0);
        assertTrue(segmenter.isWindowed());
    }
}
