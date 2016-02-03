package com.hartwig.hmftools.sullivan;

import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class OriginalFastqHeaderParserTest {

    @Test
    public void parsesCorrectly() {
        FastqHeaderParser parser = new OriginalFastqHeaderParser();
        String header = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1129:2182 1:N:0:ATCACG";
        String expectedParsed = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1129:2182";

        assertEquals(expectedParsed, parser.apply(header));
    }
}
