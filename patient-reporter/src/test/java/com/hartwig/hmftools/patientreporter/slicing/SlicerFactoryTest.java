package com.hartwig.hmftools.patientreporter.slicing;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class SlicerFactoryTest {

    private static final String BED_FILE_BASE_PATH = Resources.getResource("bed").getPath();
    private static final String UNSORTED_BED = "unsorted.bed";
    private static final String INVALID_BED = "invalid.bed";

    private static final int CPCT_GENE_PANEL_REGION_COUNT = 391;
    private static final int CPCT_GENE_PANEL_NUMBER_OF_BASES = 77840;

    private static final long GIAB_HIGH_CONFIDENCE_REGION_COUNT = 2120622L;
    private static final long GIAB_HIGH_CONFIDENCE_NUMBER_OF_BASES = 2531285550L;

    @Test
    public void canCreateCPCTGenePanelSlicer() throws IOException, EmptyFileException {
        final Slicer slicer = SlicerFactory.cpctGenePanelSlicer();
        assertEquals(CPCT_GENE_PANEL_REGION_COUNT, slicer.numberOfRegions());
        assertEquals(CPCT_GENE_PANEL_NUMBER_OF_BASES, slicer.numberOfBases());
    }

    @Test
    public void canCreateGIABHighConfidenceSlicer() throws IOException, EmptyFileException {
        final Slicer slicer = SlicerFactory.giabHighConfidenceSlicer();
        assertEquals(GIAB_HIGH_CONFIDENCE_REGION_COUNT, slicer.numberOfRegions());
        assertEquals(GIAB_HIGH_CONFIDENCE_NUMBER_OF_BASES, slicer.numberOfBases());
    }

    @Test
    public void handleUnsortedBed() throws IOException, EmptyFileException {
        String bedFile = BED_FILE_BASE_PATH + File.separator + UNSORTED_BED;
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        assertEquals(2, slicer.numberOfRegions());
    }

    @Test
    public void handleInvalidBedRegion() throws IOException, EmptyFileException {
        String bedFile = BED_FILE_BASE_PATH + File.separator + INVALID_BED;
        final Slicer slicer = SlicerFactory.fromBedFile(bedFile);
        assertEquals(2, slicer.numberOfRegions());
    }
}