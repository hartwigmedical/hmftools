package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class OutputFileUtil {

    private OutputFileUtil() {
    }

    @NotNull
    public static String generateOutputFileNameForReport(@NotNull PatientReport report) {
        SampleReport sampleReport = report.sampleReport();
        LimsCohortConfigData cohort = report.sampleReport().cohort();

        String filePrefix =
                cohort.outputFile() ? sampleReport.tumorSampleId() + "_"
                        + sampleReport.hospitalPatientId().replace(" ", "_") : sampleReport.tumorSampleId();

        String fileSuffix = report.isCorrectedReport() ? "_corrected.pdf" : ".pdf";

        String failPrefix = report instanceof QCFailReport ? "_failed" : Strings.EMPTY;

        return filePrefix + failPrefix + "_dna_analysis_report" + fileSuffix;
    }
}
