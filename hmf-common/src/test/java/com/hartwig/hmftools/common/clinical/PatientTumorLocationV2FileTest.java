package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class PatientTumorLocationV2FileTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("clinical").getPath();
    private static final String TEST_TSV = BASE_RESOURCE_DIR + File.separator + "patient_tumor_locations_v2.tsv";

    @Test
    public void canConvertDOIDs() {
        String doid1 = "doid1";
        String doid2 = "doid2";

        String doidString = PatientTumorLocationV2File.fromDOIDs(Lists.newArrayList(doid1, doid2));
        List<String> convertedDoids = PatientTumorLocationV2File.toDOIDs(doidString);

        assertEquals(2, convertedDoids.size());
        assertEquals(doid1, convertedDoids.get(0));
        assertEquals(doid2, convertedDoids.get(1));
    }

    @Test
    public void canReadFileAndConvertLines() throws IOException {
        List<PatientTumorLocationV2> patientTumorLocations = PatientTumorLocationV2File.read(TEST_TSV);

        assertEquals(3, patientTumorLocations.size());

        List<PatientTumorLocationV2> convertedTumorLocations =
                PatientTumorLocationV2File.fromLines(PatientTumorLocationV2File.toLines(patientTumorLocations));

        assertEquals(3, convertedTumorLocations.size());
        for (int i = 0; i < convertedTumorLocations.size(); i++) {
            assertEquals(patientTumorLocations.get(i), convertedTumorLocations.get(i));
        }
    }
}