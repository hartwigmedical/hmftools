package com.hartwig.hmftools.retentionchecker;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class OriginalFastqHeaderNormalizerTest {

    @Test
    public void normalizesCorrectly() {
        FastqHeaderNormalizer normalizer = new OriginalFastqHeaderNormalizer();
        String header = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1129:2182 1:N:0:ATCACG";
        String normalized = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1129:2182";

        assertEquals(normalized, normalizer.apply(header));

        String convertedHeader = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1129:2182/1";
        assertEquals(normalized, normalizer.apply(convertedHeader));
    }
}
