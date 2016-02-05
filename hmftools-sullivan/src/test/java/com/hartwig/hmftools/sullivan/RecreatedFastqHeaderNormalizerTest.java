package com.hartwig.hmftools.sullivan;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class RecreatedFastqHeaderNormalizerTest {

    @Test
    public void normalizeCorrectly() {
        FastqHeaderNormalizer normalizer = new RecreatedFastqHeaderNormalizer();
        String header = "@HISEQ_HU01:89:H7YRLADXX:1:1101:12051:6390/1";
        String expectedParsed = "@HISEQ_HU01:89:H7YRLADXX:1:1101:12051:6390";

        assertEquals(expectedParsed, normalizer.apply(header));
    }
}
