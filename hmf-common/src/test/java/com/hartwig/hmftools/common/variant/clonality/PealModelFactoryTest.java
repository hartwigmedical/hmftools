package com.hartwig.hmftools.common.variant.clonality;

import static com.hartwig.hmftools.common.variant.clonality.WeightedPloidyHistogramTest.readResource;

import static org.junit.Assert.assertEquals;

import java.util.Date;
import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Ignore;
import org.junit.Test;

public class PealModelFactoryTest {

    @Test
    public void testOffset() {
        final PeakModelFactory victim = new PeakModelFactory(10, 0.05);
        assertEquals(0.02, victim.offset(0.07), 0.01);
        assertEquals(-0.02, victim.offset(0.08), 0.01);
    }

    @Test
    public void testMaxBucket() {
        final PeakModelFactory victim = new PeakModelFactory(10, 0.05);
        victim.modelPeakHistogram(8.18, Lists.newArrayList(WeightedPloidyHistogramTest.create(8.18, 18, 55)));
    }

    @Ignore
    @Test
    public void testPeakModelling() {
        long startTime = new Date().getTime();
        WeightedPloidyHistogram victim = new WeightedPloidyHistogram(10, 0.01);
        List<ModifiableWeightedPloidy> ploidies = readResource("ploidies.tsv");

        victim.peakPloidy(10, ploidies);

        System.out.println();
        PeakModelFactory factory = new PeakModelFactory(10, 0.05);
        factory.model(ploidies);
        System.out.println(new Date().getTime() - startTime);
    }
}
