package com.hartwig.hmftools.common.purple.copynumber.sv;

import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertTrue;

import java.util.Optional;

import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.junit.Test;

public class StructuralVariantImpliedCopyNumberTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testImpliedCopyNumber() {
        final StructuralVariantPloidy left = StructuralVariantPloidyTest.create(1, Optional.of(4d), Optional.empty());
        final StructuralVariantPloidy right = StructuralVariantPloidyTest.create(-1, Optional.empty(), Optional.of(5d));
        final PurpleCopyNumber unknownCopyNumber = PurpleDatamodelTest.createCopyNumber("1", 1, 1000, 0).build();

        final PurpleCopyNumber bothKnown =
                StructuralVariantImpliedCopyNumber.impliedCopyNumber(unknownCopyNumber, Optional.of(left), Optional.of(right));
        assertTrue(bothKnown.inferred());
        assertEquals(3.5, bothKnown.averageTumorCopyNumber(), EPSILON);

        final PurpleCopyNumber leftKnown =
                StructuralVariantImpliedCopyNumber.impliedCopyNumber(unknownCopyNumber, Optional.of(left), Optional.empty());
        assertTrue(bothKnown.inferred());
        assertEquals(3, leftKnown.averageTumorCopyNumber(), EPSILON);

        final PurpleCopyNumber rightKnown =
                StructuralVariantImpliedCopyNumber.impliedCopyNumber(unknownCopyNumber, Optional.empty(), Optional.of(right));
        assertTrue(bothKnown.inferred());
        assertEquals(4, rightKnown.averageTumorCopyNumber(), EPSILON);
    }

}
