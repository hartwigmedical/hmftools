package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class PatientReporterApplicationTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();
    private static final String BED_DIRECTORY = Resources.getResource("bed").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException {
        final String bedFile = BED_DIRECTORY + File.separator + "valid.bed";
        PatientReporterApplication app = new PatientReporterApplication(RUN_DIRECTORY, bedFile, bedFile);
        app.runPatientReporter();
    }
}