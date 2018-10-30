package com.hartwig.hmftools.common.slicing;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SlicerFactoryTest {

    private static final String VALID_BED = Resources.getResource("bed/valid.bed").getPath();
    private static final String UNSORTED_BED = Resources.getResource("bed/unsorted.bed").getPath();
    private static final String INVALID_BED = Resources.getResource("bed/invalid.bed").getPath();

    @Test
    public void handleTrivialBed() throws IOException {
        final Slicer slicer = SlicerFactory.fromBedFile(VALID_BED);
        assertTrue(slicer.includes(new TestGenomePosition("1", 1)));
        assertFalse(slicer.includes(new TestGenomePosition("1", 2)));
    }

    @Test
    public void handleUnsortedBed() throws IOException {
        final Slicer slicer = SlicerFactory.fromBedFile(UNSORTED_BED);
        assertEquals(2, slicer.regions().size());
    }

    @Test
    public void handleInvalidBedRegion() throws IOException {
        final Slicer slicer = SlicerFactory.fromBedFile(INVALID_BED);
        assertEquals(2, slicer.regions().size());
    }

    private static class TestGenomePosition implements GenomePosition {

        @NotNull
        private final String chromosome;
        private final long position;

        TestGenomePosition(@NotNull final String chromosome, final long position) {
            this.chromosome = chromosome;
            this.position = position;
        }

        @NotNull
        @Override
        public String chromosome() {
            return chromosome;
        }

        @Override
        public long position() {
            return position;
        }
    }
}