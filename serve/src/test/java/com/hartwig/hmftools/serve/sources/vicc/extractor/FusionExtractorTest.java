package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.*;

import org.junit.Test;

public class FusionExtractorTest {

    @Test
    public void canExtractCodon() {
        assertEquals(Integer.valueOf(11), FusionExtractor.extractExonNumber("EXON 11 MUTATION"));
        assertEquals(Integer.valueOf(14), FusionExtractor.extractExonNumber("EXON 14 SKIPPING MUTATION"));
    }

}