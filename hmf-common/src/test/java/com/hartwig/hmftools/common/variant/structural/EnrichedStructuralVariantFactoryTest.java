package com.hartwig.hmftools.common.variant.structural;

import static com.hartwig.hmftools.common.purple.PurpleDatamodelTest.svLegPloidy;
import static com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory.adjustedCopyNumber;
import static com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory.adjustedCopyNumberChange;

import static org.junit.Assert.assertEquals;

import java.util.Optional;

import org.junit.Test;

public class EnrichedStructuralVariantFactoryTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testCopyNumberChangeWithMissingCopyNumber() {
        assertEquals(0, adjustedCopyNumberChange(svLegPloidy(1, Optional.empty(), Optional.empty(), 2).build()), EPSILON);
        assertEquals(2, adjustedCopyNumberChange(svLegPloidy(1, Optional.of(2D), Optional.empty(), 2).build()), EPSILON);
        assertEquals(2, adjustedCopyNumberChange(svLegPloidy(-1, Optional.empty(), Optional.of(2D), 2).build()), EPSILON);
    }

    @Test
    public void testMissingCopyNumber() {
        assertEquals(0, adjustedCopyNumber(svLegPloidy(1, Optional.empty(), Optional.empty(), 2).build()), EPSILON);
        assertEquals(2, adjustedCopyNumber(svLegPloidy(1, Optional.of(2D), Optional.empty(), 2).build()), EPSILON);
        assertEquals(2, adjustedCopyNumber(svLegPloidy(-1, Optional.empty(), Optional.of(2D), 2).build()), EPSILON);
    }
}
