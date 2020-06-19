package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.jetbrains.annotations.NotNull;

public final class ExtractForNumbers {

    private ExtractForNumbers() {
    }

    @NotNull
    public static String determineForNumber(@NotNull PatientReport patientReport) {
        if (patientReport instanceof QCFailReport) {
            return ((QCFailReport) patientReport).reason().forNumber();
        } else {
            assert patientReport instanceof AnalysedPatientReport;

            AnalysedPatientReport analysedPatientReport = (AnalysedPatientReport) patientReport;
            return (analysedPatientReport.hasReliablePurity() && analysedPatientReport.impliedPurity() > 0.195)
                    ? ForNumber.FOR_080.display()
                    : ForNumber.FOR_103.display();
        }
    }
}
