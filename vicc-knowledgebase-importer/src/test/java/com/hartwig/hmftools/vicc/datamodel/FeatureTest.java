package com.hartwig.hmftools.vicc.datamodel;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FeatureTest {

    @Test
    public void canDeriveProteinAnnotation() {
        assertEquals("E709K", createFeatureWithName("E709K").proteinAnnotation());
        assertEquals("E709K", createFeatureWithName("EGFR E709K").proteinAnnotation());
        assertEquals("KIT", createFeatureWithName("KIT ").proteinAnnotation());
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
