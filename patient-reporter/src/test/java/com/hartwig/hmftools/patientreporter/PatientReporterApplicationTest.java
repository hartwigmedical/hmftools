package com.hartwig.hmftools.patientreporter;

import org.junit.Test;

public class PatientReporterApplicationTest {

    @Test
    public void canRunOnRunDirectory() {
        PatientReporterApplication app = new PatientReporterApplication("RunDirKodu");
        app.printRunDirectory();
    }
}