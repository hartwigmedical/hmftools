package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class OutputFileUtil {

    private OutputFileUtil() {
    }

    @NotNull
    public static String generateOutputFileNameForReport(@NotNull PatientReport report) {
        SampleReport sampleReport = report.sampleReport();
        LimsCohort cohort = report.sampleReport().cohort();

        String filePrefix =
                cohort == LimsCohort.CORE || cohort == LimsCohort.CORELR02 || cohort == LimsCohort.CORERI02 || cohort == LimsCohort.CORELR11
                        || cohort == LimsCohort.CORESC11 || cohort == LimsCohort.COREDB ? sampleReport.tumorSampleId() + "_"
                        + sampleReport.hospitalPatientId().replace(" ", "_") : sampleReport.tumorSampleId();

        String fileSuffix = report.isCorrectedReport() ? "_corrected.pdf" : ".pdf";

        String failPrefix = report instanceof QCFailReport ? "_failed" : Strings.EMPTY;

        return filePrefix + failPrefix + "_dna_analysis_report" + fileSuffix;
    }
}
