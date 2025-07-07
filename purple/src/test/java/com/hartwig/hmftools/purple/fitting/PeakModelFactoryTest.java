package com.hartwig.hmftools.purple.fitting;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.fittingsnv.PeakModelFactory;
import com.hartwig.hmftools.purple.fittingsnv.WeightedPloidy;

import org.junit.Test;

public class PeakModelFactoryTest
{
    @Test
    public void testOffset()
    {
        PeakModelFactory peakModelFactory = new PeakModelFactory(10, 0.05);
        assertEquals(0.02, peakModelFactory.offset(0.07), 0.01);
        assertEquals(-0.02, peakModelFactory.offset(0.08), 0.01);
    }

    @Test
    public void testMaxBucket()
    {
        PeakModelFactory peakModelFactory = new PeakModelFactory(10, 0.05);
        peakModelFactory.modelPeakHistogram(8.18, Lists.newArrayList(WeightedPloidy.create(8.18, 18, 55)));
    }

    @Test
    public void testPeakLikelihood()
    {
        PeakModelFactory peakModelFactory = new PeakModelFactory(10, 0.05);
        WeightedPloidy ploidy = WeightedPloidy.create(2, 35, 50);
        assertEquals(0.06, peakModelFactory.ploidyLikelihood(1.8, ploidy), 0.001);
    }
}
