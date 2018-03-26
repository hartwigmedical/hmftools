package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ClonalityTest {

    @Test
    public void testClonality() {
        final String chromosome = "1";
        final PurityAdjuster purityAdjuster = new PurityAdjuster(Gender.FEMALE, 0.3, 1);

        assertEquals(Clonality.SUBCLONAL, Clonality.fromSample(purityAdjuster, chromosome, 3, create(4, 100)));
        assertEquals(Clonality.CLONAL, Clonality.fromSample(purityAdjuster, chromosome, 3, create(5, 100)));
        assertEquals(Clonality.CLONAL, Clonality.fromSample(purityAdjuster, chromosome, 3, create(60, 100)));
        assertEquals(Clonality.INCONSISTENT, Clonality.fromSample(purityAdjuster, chromosome, 3, create(61, 100)));
    }

    @NotNull
    private static AllelicDepth create(final int alleleReadCount, final int totalReadCount) {
        return ImmutableAllelicDepthImpl.builder().totalReadCount(totalReadCount).alleleReadCount(alleleReadCount).build();
    }
}
