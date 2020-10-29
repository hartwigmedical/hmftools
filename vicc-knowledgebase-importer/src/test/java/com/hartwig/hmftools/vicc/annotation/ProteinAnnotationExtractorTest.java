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
        assertEquals("E709K", ProteinAnnotationExtractor.proteinAnnotation(createFeatureWithName("E709K")));
        assertEquals("E709K", ProteinAnnotationExtractor.proteinAnnotation(createFeatureWithName("EGFR E709K")));
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