package com.hartwig.hmftools.retentionchecker;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class RecreatedFastqHeaderNormalizerTest {

    @Test
    public void normalizeCorrectly() {
        FastqHeaderNormalizer normalizer = new RecreatedFastqHeaderNormalizer();
        String header = "@HISEQ_HU01:89:H7YRLADXX:1:1101:12051:6390/1";
        String normalized = "@HISEQ_HU01:89:H7YRLADXX:1:1101:12051:6390";

        assertEquals(normalized, normalizer.apply(header));
    }
}
