package com.hartwig.hmftools.patientreporter.slicing;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class SlicerFactoryTest {

    private static final String BED_FILE_BASE_PATH = Resources.getResource("bed").getPath();
    private static final String FILTER_BED = "filter.bed";
    private static final String VALID_BED = "valid.bed";
    private static final String UNSORTED_BED = "unsorted.bed";
    private static final String INVALID_BED = "invalid.bed";

    @Test
    public void canCreateFilterBed() throws IOException, EmptyFileException {
        final String bedFile = BED_FILE_BASE_PATH + File.separator + FILTER_BED;
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        assertEquals(118, slicer.numberOfRegions());
        assertEquals(15593430, slicer.numberOfBases());
    }

    @Test
    public void handleTrivialBed() throws IOException, EmptyFileException {
        final String bedFile = BED_FILE_BASE_PATH + File.separator + VALID_BED;
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        assertEquals(1, slicer.numberOfRegions());
        assertEquals(1, slicer.numberOfBases());
    }

    @Test
    public void handleUnsortedBed() throws IOException, EmptyFileException {
        final String bedFile = BED_FILE_BASE_PATH + File.separator + UNSORTED_BED;
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        assertEquals(2, slicer.numberOfRegions());
    }

    @Test
    public void handleInvalidBedRegion() throws IOException, EmptyFileException {
        final String bedFile = BED_FILE_BASE_PATH + File.separator + INVALID_BED;
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        assertEquals(2, slicer.numberOfRegions());
    }
}