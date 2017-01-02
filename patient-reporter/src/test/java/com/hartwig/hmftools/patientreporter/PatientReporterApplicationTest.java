package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class PatientReporterApplicationTest {

    private static final String RUN_DIRECTORY = Resources.getResource("example").getPath();

    @Test
    public void canRunOnRunDirectory() throws IOException, HartwigException {
        PatientReporterApplication app = new PatientReporterApplication(RUN_DIRECTORY);
        app.runPatientReporter();
    }
}