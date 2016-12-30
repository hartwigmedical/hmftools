package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.junit.Test;

public class GermlineVariantFactoryTest {

    @Test
    public void extractFromNormalVCFLine() {
        String somewhatValidVCFLine = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t 10 \t 11";
        GermlineVariant variant = GermlineVariantFactory.fromVCFLine(somewhatValidVCFLine);
        assertNotNull(variant);
        assertNotNull(variant.refData());
        assertNotNull(variant.tumorData());
    }

    @Test
    public void extractFromInvalidVCFLine() {
        String invalidVCFLine = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 ";
        assertNull(GermlineVariantFactory.fromVCFLine(invalidVCFLine));
    }
}