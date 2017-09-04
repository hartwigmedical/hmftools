package com.hartwig.hmftools.common.ecrf.doid;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.EmptyFileException;

import org.junit.Test;

public class DoidMappingTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("ecrf").getPath();
    private static final String TEST_FILE =
            BASE_RESOURCE_DIR + File.separator + "doid" + File.separator + "tumor_location_doid_mapping.csv";

    @Test
    public void canLoadFromCsv() throws IOException, EmptyFileException {

        final Map<String, Set<String>> doidMapping = TumorLocationDoidMapping.readMappingFromCSV(TEST_FILE);
        assertEquals(72, doidMapping.size());
        assertEquals(2, doidMapping.get("Breast Cancer: ER-negative/HER2-positive").size());
        assertTrue(doidMapping.get("Breast Cancer: ER-negative/HER2-positive").contains("DOID:0060076"));
        assertTrue(doidMapping.get("Breast Cancer: ER-negative/HER2-positive").contains("DOID:0060079"));
        assertEquals(1, doidMapping.get("Eye cancer").size());
        assertTrue(doidMapping.get("Eye cancer").contains("DOID:2174"));
        assertEquals(0, doidMapping.get("Breast Cancer: subtype unknown").size());
    }
}
