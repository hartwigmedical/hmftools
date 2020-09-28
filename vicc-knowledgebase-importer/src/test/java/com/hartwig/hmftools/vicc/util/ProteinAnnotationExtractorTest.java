package com.hartwig.hmftools.vicc.util;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ProteinAnnotationExtractorTest {

    @Test
    public void canExtractProteinAnnotationFromFeature() {
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation(createFeatureWithName("E709K")));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation(createFeatureWithName("EGFR E709K")));
        assertEquals("KIT", ProteinAnnotationExtractor.toProteinAnnotation(createFeatureWithName("KIT ")));
    }

    @Test
    public void canConvertFeatureNameToProteinAnnotation() {
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR E709K "));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR:E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR:p.E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("EGFR p.E709K"));
        assertEquals("E709K", ProteinAnnotationExtractor.toProteinAnnotation("E709K (c.2100A>c)"));

        assertEquals("G778_P780dup", ProteinAnnotationExtractor.toProteinAnnotation("G778_P780DUP"));
        assertEquals("V560del", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560DEL"));
        assertEquals("V560fs", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560FS"));
        assertEquals("V560fs", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560FS*"));
        assertEquals("V560insAYVM", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560INSAYVM"));
        assertEquals("V560insINS", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560INSINS"));
        assertEquals("V560delinsDEL", ProteinAnnotationExtractor.toProteinAnnotation("KIT:p.V560DELINSDEL"));
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