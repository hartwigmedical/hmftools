package com.hartwig.hmftools.iclusion.classification;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class IclusionProteinAnnotationExtractorTest {

    @Test
    public void canGenerateProteinAnnotation() {
        IclusionProteinAnnotationExtractor extractor = new IclusionProteinAnnotationExtractor();
        assertEquals("W288fs", extractor.apply("W288FS"));
        assertEquals("W288IFS", extractor.apply("W288IFS"));

        assertEquals("W288delinsDEL", extractor.apply("W288DELINSDEL"));
    }
}