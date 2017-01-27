package com.hartwig.hmftools.retentionchecker;

import htsjdk.samtools.fastq.FastqRecord;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class FastqHeaderTest {

    @Test
    public void createCorrectlyFromFastqRecord() {
        final String headerString = "@HISEQ_HU01:89:H7YRLADXX:1:1101:1129:2182";
        final FastqRecord record = new FastqRecord(headerString, "", null, "");

        final FastqHeader header = FastqHeader.parseFromFastqRecord(record, new DummyNormalizer());

        final FastqHeader expectedHeader = new FastqHeader("@HISEQ_HU01", 89, "H7YRLADXX", 1, 1101, 1129, 2182);

        assertTrue(expectedHeader.equals(header));
    }

    private static final class DummyNormalizer implements FastqHeaderNormalizer {
        public String apply(String s) {
            return s;
        }
    }
}
