package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.nio.file.NoSuchFileException;

import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class PatientReporterApplicationTest {

    @Test(expected = NoSuchFileException.class)
    public void canRunOnRunDirectory() throws IOException, HartwigException {
        PatientReporterApplication app = new PatientReporterApplication("RunDirKodu");
        app.runPatientReporter();
    }
}