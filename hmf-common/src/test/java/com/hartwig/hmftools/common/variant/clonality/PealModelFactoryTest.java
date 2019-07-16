package com.hartwig.hmftools.common.variant.clonality;

import static com.hartwig.hmftools.common.variant.clonality.WeightedPloidyHistogramTest.readResource;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Date;
import java.util.List;

import org.junit.Ignore;
import org.junit.Test;

public class PealModelFactoryTest {

    @Test
    public void testOffset() {
        final PeakModelFactory victim = new PeakModelFactory(10, 0.05);
        assertEquals(0.02, victim.offset(0.07), 0.01);
        assertEquals(-0.02, victim.offset(0.08), 0.01);
    }

    @Ignore
    public void testPeakModelling() throws IOException {
        long startTime = new Date().getTime();
        final WeightedPloidyHistogram victim = new WeightedPloidyHistogram(10, 0.01);
        List<ModifiableWeightedPloidy> ploidies = readResource("ploidies.tsv");
        double[] histogram = victim.histogram(ploidies);

        double peakPloidy = victim.peakPloidy(10, histogram);

        System.out.println();
        PeakModelFactory factory = new PeakModelFactory(10, 0.05);
        List<PeakModel> result = factory.model(ploidies);
        System.out.println(new Date().getTime() - startTime);
    }

}
