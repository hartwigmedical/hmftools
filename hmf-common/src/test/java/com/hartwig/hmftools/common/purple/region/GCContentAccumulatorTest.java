package com.hartwig.hmftools.common.purple.region;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.copynumber.freec.FreecGCContent;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;

import org.junit.Before;
import org.junit.Test;

public class GCContentAccumulatorTest {

    private static final double EPSILON = 1e-10;

    private ObservedRegionFactory.GCContentAccumulator victim;

    @Before
    public void setup() {
        victim = new ObservedRegionFactory.GCContentAccumulator();
        assertGCContent(0, 0, 0);
    }

    @Test
    public void testStandardBehaviour() {
        victim.accept(create(0.8, 0.6, 1));
        assertGCContent(0.8, 0.6, 1);

        victim.accept(create(0.9, 0.9, 0.8));
        assertGCContent(0.86, 0.75, 0.9);
    }

    @Test
    public void testNoGcContent() {
        victim.accept(create(-1, 0, 0));
        assertGCContent(0, 0, 0);

        victim.accept(create(1, 1, 1));
        assertGCContent(1, 0.5, 0.5);
    }

    private void assertGCContent(double gcContent, double nonNPercentage, double mappablePercentage) {
        assertEquals(gcContent, victim.getAverageGCContent(), EPSILON);
        assertEquals(nonNPercentage, victim.getAverageNonNPercentage(), EPSILON);
        assertEquals(mappablePercentage, victim.getSumMappablePercentage(), EPSILON);
    }

    private static FreecGCContent create(double gcContent, double nonNPercentage, double mappablePercentage) {
        return PurpleDatamodelTest.createGCContent("1", 1, gcContent, nonNPercentage, mappablePercentage).build();
    }
}
