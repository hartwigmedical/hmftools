package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ClonalityTest {

    @Test
    public void testClonality() {
        assertEquals(Clonality.SUBCLONAL, Clonality.fromSample(3, 0.3, create(4, 100)));
        assertEquals(Clonality.CLONAL, Clonality.fromSample(3, 0.3, create(5, 100)));
        assertEquals(Clonality.CLONAL, Clonality.fromSample(3, 0.3, create(60, 100)));
        assertEquals(Clonality.INCONSISTENT, Clonality.fromSample(3, 0.3, create(61, 100)));
    }


    private static AllelicDepth create(final int alleleReadCount, final int totalReadCount) {
        return new AllelicDepth() {
            @Override
            public int totalReadCount() {
                return totalReadCount;
            }

            @Override
            public int alleleReadCount() {
                return alleleReadCount;
            }
        };
    }

}
