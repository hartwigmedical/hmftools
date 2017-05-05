package com.hartwig.hmftools.patientreporter.slicing;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import static org.junit.Assert.*;

public class HMFSlicingAnnotationFactoryTest {

    @Test
    public void canCreateFromGenomeRegion() {
        assertNull(HMFSlicingAnnotationFactory.fromGenomeRegion(toRegion(null)));
        assertNull(HMFSlicingAnnotationFactory.fromGenomeRegion(toRegion("korneel")));
        assertNull(HMFSlicingAnnotationFactory.fromGenomeRegion(toRegion("hello world")));

        final HMFSlicingAnnotation annotation = HMFSlicingAnnotationFactory.fromGenomeRegion(
                toRegion("ENST00000361445.4 (MTOR)"));
        assertNotNull(annotation);
        assertEquals("ENST00000361445", annotation.transcriptID());
        assertEquals(4, annotation.transcriptVersion());
        assertEquals("ENST00000361445.4", annotation.transcript());
        assertEquals("MTOR", annotation.gene());
    }

    @NotNull
    private static GenomeRegion toRegion(@Nullable final String annotation) {
        return ImmutableBEDGenomeRegion.of(Strings.EMPTY, 0L, 0L, annotation);
    }
}