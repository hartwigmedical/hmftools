package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class CkbProteinAnnotationExtractorTest {

    @Test
    public void canExtractProteinAnnotationFromEvent() {
        CkbProteinAnnotationExtractor extractor = new CkbProteinAnnotationExtractor();

        assertEquals("E709K", extractor.apply("E709K"));
        assertEquals("E709fs", extractor.apply("E709fs*2"));
    }
}