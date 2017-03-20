package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class VariantConsequenceTest {

    @Test
    public void canDetectSubTypes() {
        assertTrue(VariantConsequence.FRAMESHIFT_VARIANT.isParentTypeOf("frameshift_variant"));
        assertTrue(VariantConsequence.FRAMESHIFT_VARIANT.isParentTypeOf("plus_2_frameshift_variant"));
        assertFalse(VariantConsequence.FRAMESHIFT_VARIANT.isParentTypeOf("whatever"));
    }
}