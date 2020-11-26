package com.hartwig.hmftools.iclusion.classification;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ProteinAnnotationExtractorTest {

    @Test
    public void canGenerateProteinAnnotation() {
        ProteinAnnotationExtractor extractor = new ProteinAnnotationExtractor();
        assertEquals("W288fs", extractor.apply("W288FS"));
        assertEquals("W288IFS", extractor.apply("W288IFS"));

        assertEquals("W288delinsDEL", extractor.apply("W288DELINSDEL"));
    }
}