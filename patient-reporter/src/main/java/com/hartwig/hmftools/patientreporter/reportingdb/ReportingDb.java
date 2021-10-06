package com.hartwig.hmftools.patientreporter.reportingdb;

import static java.lang.String.format;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Map;

import com.google.gson.GsonBuilder;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.reportingdb.ReportingDatabase;
import com.hartwig.hmftools.common.reportingdb.ReportingEntry;
import com.hartwig.hmftools.patientreporter.SampleReport;
import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.algo.GenomicAnalysis;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReportingDb {

    private static final Logger LOGGER = LogManager.getLogger(ReportingDb.class);

    private static final String NA_STRING = "N/A";

    public ReportingDb() {
    }

    public void appendAnalysedReport(@NotNull AnalysedPatientReport report, @NotNull String outputDirectory) throws IOException {
        if (shouldBeAddedToReportingDb(report.sampleReport())) {
            String sampleId = report.sampleReport().tumorSampleId();
            LimsCohortConfig cohort = report.sampleReport().cohort();

            if (report.sampleReport().cohort().reportConclusion() && report.clinicalSummary().isEmpty()) {
                LOGGER.warn("Skipping addition to reporting db, missing summary for sample '{}'!", sampleId);
            } else {
                String tumorBarcode = report.sampleReport().tumorSampleBarcode();
                String reportDate = ReportResources.REPORT_DATE;
                GenomicAnalysis analysis = report.genomicAnalysis();

                String purity = new DecimalFormat("0.00").format(analysis.impliedPurity());
                boolean hasReliableQuality = analysis.hasReliableQuality();
                boolean hasReliablePurity = analysis.hasReliablePurity();

                String reportType;
                if (hasReliablePurity && analysis.impliedPurity() > ReportResources.PURITY_CUTOFF) {
                    reportType = "dna_analysis_report";
                } else {
                    reportType = "dna_analysis_report_insufficient_tcp";
                }

                if (report.isCorrectedReport()) {
                    reportType = reportType + "_corrected";
                }

                writeApiUpdateJson(outputDirectory,
                        tumorBarcode,
                        sampleId,
                        cohort,
                        reportType,
                        reportDate,
                        purity,
                        hasReliableQuality,
                        hasReliablePurity);
            }
        }
    }

    private void writeApiUpdateJson(final String outputDirectory, final String tumorBarcode, final String sampleId,
            final LimsCohortConfig cohort, final String reportType, final String reportDate, final String purity,
            final Boolean hasReliableQuality, final Boolean hasReliablePurity) throws IOException {
        File outputFile = new File(outputDirectory, format("%s_%s_%s_api-update.json", sampleId, tumorBarcode, reportType));
        Map<String, Object> payload = new HashMap<>();
        payload.put("barcode", tumorBarcode);
        payload.put("report_type", reportType);
        payload.put("report_date", reportDate);
        payload.put("purity", purity.equals(NA_STRING) ? purity : Float.parseFloat(purity));
        payload.put("cohort", cohort.cohortId());
        payload.put("has_reliable_quality", hasReliableQuality != null ? hasReliableQuality : NA_STRING);
        payload.put("has_reliable_purity", hasReliablePurity != null ? hasReliablePurity : NA_STRING);

        appendToFile(outputFile.getAbsolutePath(),
                new GsonBuilder().serializeNulls()
                        .serializeSpecialFloatingPointValues()
                        .setPrettyPrinting()
                        .disableHtmlEscaping()
                        .create()
                        .toJson(payload));
    }

    public void appendQCFailReport(@NotNull QCFailReport report, @NotNull String outputDirectory) throws IOException {
        if (shouldBeAddedToReportingDb(report.sampleReport())) {
            String sampleId = report.sampleReport().tumorSampleId();
            LimsCohortConfig cohort = report.sampleReport().cohort();
            String tumorBarcode = report.sampleReport().tumorSampleBarcode();
            String reportDate = ReportResources.REPORT_DATE;

            String reportType = report.isCorrectedReport() ? report.reason().identifier() + "_corrected" : report.reason().identifier();

            writeApiUpdateJson(outputDirectory, tumorBarcode, sampleId, cohort, reportType, reportDate, NA_STRING, null, null);

        }
    }

    private static boolean shouldBeAddedToReportingDb(@NotNull SampleReport report) {
        String sampleId = report.tumorSampleId();
        if (sampleId.startsWith("COLO")) {
            LOGGER.info("Sample '{}' filtered for reporting db because it appears to be belong to COLO test samples", sampleId);
            return false;
        } else if (report.cohort().cohortId().isEmpty()) {
            LOGGER.info("Sample '{}' filtered for reporting db since it does not belong to a cohort and likely a test sample", sampleId);
            return false;
        }
        return true;
    }

    private static void appendToFile(@NotNull String reportingDbTsv, @NotNull String stringToAppend) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(reportingDbTsv, true));
        writer.write(stringToAppend);
        writer.close();
    }
}
