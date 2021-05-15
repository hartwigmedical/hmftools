package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class ViccProteinAnnotationExtractorTest {

    @Test
    public void canExtractProteinAnnotationFromEvent() {
        ViccProteinAnnotationExtractor extractor = new ViccProteinAnnotationExtractor();
        assertEquals("E709K", extractor.apply("E709K"));
        assertEquals("E709K", extractor.apply("EGFR E709K "));
        assertEquals("E709K", extractor.apply("EGFR:E709K"));
        assertEquals("E709K", extractor.apply("EGFR:p.E709K"));
        assertEquals("E709K", extractor.apply("EGFR p.E709K"));
        assertEquals("E709K", extractor.apply("E709K (c.2100A>c)"));

        assertEquals("G778_P780dup", extractor.apply("G778_P780DUP"));
        assertEquals("V560del", extractor.apply("KIT:p.V560DEL"));
        assertEquals("V560fs", extractor.apply("KIT:p.V560FS"));
        assertEquals("V560fs", extractor.apply("KIT:p.V560FS*"));
        assertEquals("V560insAYVM", extractor.apply("KIT:p.V560INSAYVM"));
        assertEquals("V560insINS", extractor.apply("KIT:p.V560INSINS"));
        assertEquals("V560delinsDEL", extractor.apply("KIT:p.V560DELINSDEL"));

        assertEquals("S978fs", extractor.apply("S978FS*4"));
        assertEquals("S978fs", extractor.apply("S978fs*123"));
    }
}