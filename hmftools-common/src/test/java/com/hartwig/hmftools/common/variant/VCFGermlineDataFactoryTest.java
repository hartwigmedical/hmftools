package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.common.variant.VCFGermlineDataFactory;

import org.junit.Test;

public class VCFGermlineDataFactoryTest {

    @Test
    public void extractFromNormalVCF() {
        String somewhatValidVCF = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 \t 10 \t 11";
        assertNotNull(VCFGermlineDataFactory.fromVCFLine(somewhatValidVCF));
    }

    @Test
    public void extractFromInvalidVCF() {
        String somewhatValidVCF = "1 \t 2 \t 3 \t 4 \t 5 \t 6 \t 7 \t 8 \t 9 ";
        assertNull(VCFGermlineDataFactory.fromVCFLine(somewhatValidVCF));
    }
}