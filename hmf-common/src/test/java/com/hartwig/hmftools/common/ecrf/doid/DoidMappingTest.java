package com.hartwig.hmftools.common.ecrf.doid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class DoidMappingTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("ecrf").getPath();
    private static final String TEST_FILE =
            BASE_RESOURCE_DIR + File.separator + "doid" + File.separator + "tumor_location_doid_mapping.csv";

    @Test
    public void canLoadFromCsv() throws IOException, EmptyFileException {
        final TumorLocationDoidMapping doidMapping = TumorLocationDoidMapping.fromCSV(TEST_FILE);
        assertEquals(2, doidMapping.doidsForTumorType("Breast Cancer: ER-negative/HER2-positive").size());
        assertTrue(doidMapping.doidsForTumorType("Breast Cancer: ER-negative/HER2-positive").contains("DOID:0060076"));
        assertTrue(doidMapping.doidsForTumorType("Breast Cancer: ER-negative/HER2-positive").contains("DOID:0060079"));
        assertEquals(1, doidMapping.doidsForTumorType("Eye cancer").size());
        assertTrue(doidMapping.doidsForTumorType("Eye cancer").contains("DOID:2174"));
        assertEquals(0, doidMapping.doidsForTumorType("Breast Cancer: subtype unknown").size());
    }
}
