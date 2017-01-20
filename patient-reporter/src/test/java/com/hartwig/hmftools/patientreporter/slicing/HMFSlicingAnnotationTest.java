package com.hartwig.hmftools.patientreporter.slicing;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class HMFSlicingAnnotationTest {

    @Test
    public void canCreateFromGenomeRegion() {
        assertNull(HMFSlicingAnnotation.fromGenomeRegion(toRegion(null)));
        assertNull(HMFSlicingAnnotation.fromGenomeRegion(toRegion("korneel")));
        assertNull(HMFSlicingAnnotation.fromGenomeRegion(toRegion("hello world")));

        final HMFSlicingAnnotation annotation = HMFSlicingAnnotation.fromGenomeRegion(
                toRegion("ENST00000361445.4 (MTOR)"));
        assertNotNull(annotation);
        assertEquals("ENST00000361445", annotation.transcriptID());
        assertEquals(4, annotation.transcriptVersion());
        assertEquals("ENST00000361445.4", annotation.transcript());
        assertEquals("MTOR", annotation.gene());
    }

    @NotNull
    private static GenomeRegion toRegion(@Nullable final String annotation) {
        return new GenomeRegion(Strings.EMPTY, 0L, 0L, annotation);
    }
}