package com.hartwig.hmftools.vicc.util;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ProteinAnnotationExtractorTest {

    @Test
    public void canConvertNameToProteinAnnotation() {
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR E709K "));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR:E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR:p.E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR p.E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("E709K (c.2100A>c)"));

        assertEquals("V560del", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560DEL"));
        assertEquals("V560fs", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560FS"));
        assertEquals("V560insAYVM", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560INSAYVM"));
        assertEquals("V560insINS", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560INSINS"));
        assertEquals("V560delinsDEL", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560DELINSDEL"));
    }
}