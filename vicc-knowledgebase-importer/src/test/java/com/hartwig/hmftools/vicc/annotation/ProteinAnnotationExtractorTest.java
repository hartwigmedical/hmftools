package com.hartwig.hmftools.vicc.annotation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ProteinAnnotationExtractorTest {

    @Test
    public void canExtractProteinAnnotationFromFeature() {
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation(createFeatureWithName("E709K")));
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation(createFeatureWithName("EGFR E709K")));
        assertEquals("KIT", ProteinAnnotationExtractor.extractProteinAnnotation(createFeatureWithName("KIT ")));
    }

    @Test
    public void canConvertFeatureNameToProteinAnnotation() {
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation("E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation("EGFR E709K "));
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation("EGFR:E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation("EGFR:p.E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation("EGFR p.E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.extractProteinAnnotation("E709K (c.2100A>c)"));

        assertEquals("G778_P780dup", ProteinAnnotationExtractor.extractProteinAnnotation("G778_P780DUP"));
        assertEquals("V560del", ProteinAnnotationExtractor.extractProteinAnnotation("KIT:p.V560DEL"));
        assertEquals("V560fs", ProteinAnnotationExtractor.extractProteinAnnotation("KIT:p.V560FS"));
        assertEquals("V560fs", ProteinAnnotationExtractor.extractProteinAnnotation("KIT:p.V560FS*"));
        assertEquals("V560insAYVM", ProteinAnnotationExtractor.extractProteinAnnotation("KIT:p.V560INSAYVM"));
        assertEquals("V560insINS", ProteinAnnotationExtractor.extractProteinAnnotation("KIT:p.V560INSINS"));
        assertEquals("V560delinsDEL", ProteinAnnotationExtractor.extractProteinAnnotation("KIT:p.V560DELINSDEL"));
    }

    @NotNull
    private static Feature createFeatureWithName(@NotNull String name) {
        return ImmutableFeature.builder()
                .name(name)
                .biomarkerType(Strings.EMPTY)
                .referenceName(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .start(Strings.EMPTY)
                .end(Strings.EMPTY)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .build();
    }
}