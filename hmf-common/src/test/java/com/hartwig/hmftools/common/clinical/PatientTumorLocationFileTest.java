package com.hartwig.hmftools.common.clinical;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PatientTumorLocationFileTest {

    private static final String BASE_RESOURCE_DIR = Resources.getResource("clinical").getPath();
    private static final String TEST_CSV = BASE_RESOURCE_DIR + File.separator + "patient_tumor_locations.csv";
    private static final String TEST_TSV = BASE_RESOURCE_DIR + File.separator + "patient_tumor_locations.tsv";

    private static final PatientTumorLocation PATIENT1 = ImmutablePatientTumorLocation.builder()
            .patientIdentifier("patient1")
            .primaryTumorLocation("Colon/Rectum")
            .cancerSubtype(Strings.EMPTY)
            .build();

    private static final PatientTumorLocation PATIENT2 = ImmutablePatientTumorLocation.builder()
            .patientIdentifier("patient2")
            .primaryTumorLocation("Lung")
            .cancerSubtype("Non-Small Cell")
            .build();

    @Test
    public void canReadInputFiles() throws IOException {
        List<PatientTumorLocation> tumorLocationsFromCsv = PatientTumorLocationFile.readRecordsCSV(TEST_CSV);
        assertTumorLocations(tumorLocationsFromCsv);

        List<PatientTumorLocation> tumorLocationsFromTsv = PatientTumorLocationFile.readRecordsTSV(TEST_TSV);
        assertTumorLocations(tumorLocationsFromTsv);
    }

    private static void assertTumorLocations(@NotNull List<PatientTumorLocation> patientTumorLocations) {
        assertEquals(2, patientTumorLocations.size());

        PatientTumorLocation patient1 = patientTumorLocations.get(0);
        assertEquals(PATIENT1.patientIdentifier(), patient1.patientIdentifier());
        assertEquals(PATIENT1.primaryTumorLocation(), patient1.primaryTumorLocation());
        assertEquals(PATIENT1.cancerSubtype(), patient1.cancerSubtype());

        PatientTumorLocation patient2 = patientTumorLocations.get(1);
        assertEquals(PATIENT2.patientIdentifier(), patient2.patientIdentifier());
        assertEquals(PATIENT2.primaryTumorLocation(), patient2.primaryTumorLocation());
        assertEquals(PATIENT2.cancerSubtype(), patient2.cancerSubtype());
    }

}